#!/usr/bin/env python
#-*- coding: utf-8 -*-

###########################################################################
##                                                                       ##
## Copyrights Etienne Chové <chove@crans.org> 2009                       ##
## Copyrights Frédéric Rodrigo 2011-2015                                 ##
##                                                                       ##
## This program is free software: you can redistribute it and/or modify  ##
## it under the terms of the GNU General Public License as published by  ##
## the Free Software Foundation, either version 3 of the License, or     ##
## (at your option) any later version.                                   ##
##                                                                       ##
## This program is distributed in the hope that it will be useful,       ##
## but WITHOUT ANY WARRANTY; without even the implied warranty of        ##
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         ##
## GNU General Public License for more details.                          ##
##                                                                       ##
## You should have received a copy of the GNU General Public License     ##
## along with this program.  If not, see <http://www.gnu.org/licenses/>. ##
##                                                                       ##
###########################################################################

from Analyser_Osmosis import Analyser_Osmosis

import os
import math
import numpy
import timeit
import pickle
from shapely.geometry.linestring    import LineString


sql10 = """
CREATE TEMP TABLE {0}buildings AS
SELECT
    ways.id,
    ways.linestring,
    ways.tags->'building' AS building,
    NOT ways.tags?'wall' OR ways.tags->'wall' != 'no' AS wall,
    array_length(ways.nodes, 1) AS nodes_length,
    ST_MakePolygon(ways.linestring) AS polygon
FROM
    {0}ways AS ways
    LEFT JOIN relation_members ON
        relation_members.member_id = ways.id AND
        relation_members.member_type = 'W'
WHERE
    relation_members.member_id IS NULL AND
    (ways.tags?'building' AND ways.tags->'building' != 'no') AND
    NOT ways.tags?'layer' AND
    is_polygon AND
    ST_IsValid(ways.linestring) = 't' AND
    ST_IsSimple(ways.linestring) = 't'
"""

sql11 = """
CREATE INDEX {0}buildings_polygon_idx ON {0}buildings USING gist(polygon);
CREATE INDEX {0}buildings_wall_idx ON {0}buildings(wall);
"""

sql20 = """
CREATE TEMP TABLE {0}bnodes AS
SELECT
    id,
    ST_PointN(linestring, generate_series(1, ST_NPoints(linestring))) AS geom
FROM
    {0}buildings
WHERE
    wall
"""

sql21 = """
CREATE INDEX {0}bnodes_geom ON {0}bnodes USING GIST(geom);
"""

sql30 = """
CREATE TABLE intersection_{0}_{1} AS
SELECT
    b1.id AS id1,
    b2.id AS id2,
    ST_AsText(ST_Centroid(ST_Intersection(b1.polygon, b2.polygon))),
    ST_Area(ST_Intersection(b1.polygon, b2.polygon)) AS intersectionArea,
    least(ST_Area(b1.polygon), ST_Area(b2.polygon))*0.10 AS threshold
FROM
    {0}buildings AS b1,
    {1}buildings AS b2
WHERE
    b1.id > b2.id AND
    b1.wall AND
    b2.wall AND
    b1.polygon && b2.polygon AND
    ST_Area(ST_Intersection(b1.polygon, b2.polygon)) <> 0
"""

sql31 = """
SELECT
    *
FROM
    intersection_{0}_{1}
"""

sql40 = """
SELECT
    id,
    ST_AsText(ST_Centroid(polygon))
FROM
    {0}buildings
WHERE
    wall AND
    ST_Area(polygon) < 0.05e-10
"""

sql50 = """
SELECT
    DISTINCT ON (bnodes.id, bnodes.geom)
    buildings.id,
    bnodes.id,
    ST_AsText(bnodes.geom)
FROM
    {0}bnodes AS bnodes
    JOIN {1}buildings AS buildings ON
        buildings.id != bnodes.id AND
        buildings.wall AND
        ST_DWithin(buildings.linestring, bnodes.geom, 1e-7) AND
        ST_Distance(buildings.linestring, bnodes.geom) > 0
ORDER BY
    bnodes.id,
    bnodes.geom
"""

sql60 = """
SELECT
    ST_AsText(ST_Centroid(geom)),
    ST_Area(geom)
FROM
    (
    SELECT
        (ST_Dump(poly)).geom AS geom
    FROM
        (
        SELECT
            ST_Union(ST_Buffer(ways.linestring,5e-3,'quad_segs=2')) AS poly
        FROM
            intersection_{0}_{1}
            JOIN ways ON
                ways.id = id1
        WHERE
            intersectionArea > threshold
        ) AS building
    ) AS buffer
WHERE
    ST_Area(geom) > 5e-4
"""

sql70 = """
SELECT
   DISTINCT ON (b2.id)
   b2.id,
   ST_AsText(way_locate(b2.linestring))
FROM
   {0}buildings AS b1,
   {1}buildings AS b2
WHERE
   b2.nodes_length = 4 AND
   b2.id != b1.id AND
   b1.building = b2.building AND
   b1.wall = b2.wall AND
   ST_Intersects(b1.polygon, b2.polygon)
ORDER BY
   b2.id
"""

sql80 = """
SELECT
   b1.id,
   b2.id,
   ST_AsText(way_locate(b2.linestring)),
   ST_AsText(ST_Transform(b1.polygon, {2})),
   ST_AsText(ST_Transform(b2.polygon, {2}))
FROM
   {0}buildings AS b1,
   {1}buildings AS b2
WHERE
   b1.id < b2.id AND
   b1.building = b2.building AND
   b1.wall = b2.wall AND
   ST_Intersects(b1.polygon, b2.polygon)
"""

class Analyser_Osmosis_Building_Overlaps(Analyser_Osmosis):

    def __init__(self, config, logger = None):
        Analyser_Osmosis.__init__(self, config, logger)
        self.FR = config.options and ("country" in config.options and config.options["country"] == "FR" or "test" in config.options)
        self.classs_change[1] = {"item":"0", "level": 3, "tag": ["building", "geom", "fix:chair"], "desc": T_(u"Building intersection") }
        self.classs_change[2] = {"item":"0", "level": 2, "tag": ["building", "geom", "fix:chair"], "desc": T_(u"Large building intersection") }
        self.classs_change[3] = {"item":"0", "level": 3, "tag": ["building", "geom", "fix:chair"], "desc": T_(u"Building too small") }
        self.classs_change[4] = {"item":"0", "level": 3, "tag": ["building", "geom", "fix:chair"], "desc": T_(u"Gap between buildings") }
        self.classs_change[5] = {"item":"0", "level": 1, "tag": ["building", "fix:chair"], "desc": T_(u"Large building intersection cluster") }
        if self.FR:
            self.classs_change[6] = {"item":"1", "level": 3, "tag": ["building", "geom", "fix:chair"], "desc": T_(u"Building in parts") }
            self.classs_change[7] = {"item":"1", "level": 3, "tag": ["building", "geom", "fix:chair"], "desc": T_(u"Building in parts") }
        self.callback30 = lambda res: {"class":2 if res[3]>res[4] else 1, "data":[self.way, self.way, self.positionAsText]}
        self.callback40 = lambda res: {"class":3, "data":[self.way, self.positionAsText]}
        self.callback50 = lambda res: {"class":4, "data":[self.way, self.way, self.positionAsText]}
        self.callback60 = lambda res: {"class":5, "data":[self.positionAsText]}
        if self.FR:
            self.callback70 = lambda res: {"class":6, "data":[self.way, self.positionAsText]}
            self.callback80 = lambda res: {"class":7, "data":[self.way, self.way, self.positionAsText]} \
                    if self.segmentedDetector.check(res[3],res[4]) else None
            self.segmentedDetector = SegmentedDetector()

    def analyser_osmosis_all(self):
        self.run(sql10.format(""))
        self.run(sql11.format(""))
        self.run(sql20.format(""))
        self.run(sql21.format(""))
        self.run(sql30.format("", ""))
        self.run(sql31.format("", ""), self.callback30)
        self.run(sql40.format(""), self.callback40)
        self.run(sql50.format("", ""), self.callback50)
        self.run(sql60.format("", ""), self.callback60)
        if self.FR:
            self.run(sql70.format("", ""), self.callback70)
            t = Timer("sql80")
            self.run(sql80.format("", "", self.config.options.get("proj")), self.callback80)
            t.prnt()

    def analyser_osmosis_touched(self):
        self.run(sql10.format(""))
        self.run(sql11.format(""))
        self.run(sql20.format(""))
        self.run(sql21.format(""))
        self.run(sql10.format("touched_"))
        self.run(sql11.format("touched_"))
        self.run(sql20.format("touched_"))
        self.run(sql21.format("touched_"))
        dup = set()
        self.run(sql30.format("touched_", ""))
        self.run(sql30.format("", "touched_"))
        self.run(sql30.format("touched_", "touched_"))
        self.run(sql31.format("touched_", ""), lambda res: dup.add(res[0]) or self.callback30(res))
        self.run(sql31.format("", "touched_"), lambda res: res[0] in dup or dup.add(res[0]) or self.callback30(res))
        self.run(sql31.format("touched_", "touched_"), lambda res: res[0] in dup or dup.add(res[0]) or self.callback30(res))
        self.run(sql40.format("touched_"), self.callback40)
        dup = set()
        self.run(sql50.format("touched_", ""), lambda res: dup.add(res[0]) or self.callback50(res))
        self.run(sql50.format("", "touched_"), lambda res: res[0] in dup or dup.add(res[0]) or self.callback50(res))
        self.run(sql50.format("touched_", "touched_"), lambda res: res[0] in dup or dup.add(res[0]) or self.callback50(res))
        #self.run(sql60.format("", ""), self.callback60) Can be done in diff mode without runing a full sql30
        if self.FR:
            self.run(sql70.format("touched_", ""), self.callback70)
            self.run(sql70.format("", "touched_"), self.callback70)
            self.run(sql70.format("touched_", "touched_"), self.callback70)
            self.run(sql80.format("touched_", "", self.config.options.get("proj")), self.callback80)
            self.run(sql80.format("", "touched_", self.config.options.get("proj")), self.callback80)
            self.run(sql80.format("touched_", "touched_", self.config.options.get("proj")), self.callback80)


class SegmentedDetector(object):
    """Detector for segmented building due to French cadastre artifacts."""
    def __init__(self):
        self.classifier = SegmentedDetector.load_classifier()
    def check(self, p1, p2):
        #Quick parse of the POLYGON((x1 y1, x2 y2, ...))" WKT string:
        p1 = map(lambda e:map(float,e.split(" ")), p1[9:-2].split(","))
        p2 = map(lambda e:map(float,e.split(" ")), p2[9:-2].split(","))
        vector1 = SegmentedDetector.get_segmented_analysis_vector(p1, p2)
        vector2 = SegmentedDetector.get_segmented_analysis_vector(p2, p1)
        result =  ((vector1 != None) and (self.classifier.predict(vector1) == 1).all()) or \
                ((vector2 != None) and (self.classifier.predict(vector2) == 1).all())
        return result

    @staticmethod
    def load_classifier():
        #FIXME: use config instead of hardcoded path:
        segmented_data_dir="/home/osm/cadastre-housenumber/segmented_building_data"
        #os.system("cd " + segmented_data_dir  +"; make -s")
        classifier = pickle.load(open(os.path.join(segmented_data_dir, "classifier.pickle")))
        return classifier

    @staticmethod
    def get_segmented_analysis_vector_from_polygons(p1, p2):
        assert(len(p1.interiors) == 0)
        assert(len(p2.interiors) == 0)
        return SegmentedDetector.get_segmented_analysis_vector(p1.exterior.coords, p2.exterior.coords)

    @staticmethod
    def get_segmented_analysis_vector(way1, way2):
        result = None
        if (way1[-1] == way1[0]) and (way2[-1] == way2[0]):
            external1, common, external2 = SegmentedDetector.get_external1_common_external2_ways(way1, way2)
            if len(common)>1:
                assert(external1[-1] == common[0])
                assert(external2[-1] == common[0])
                assert(common[-1] == external1[0])
                assert(common[-1] == external2[0])

                #        a-----------b-------------c
                #        |            \            |
                #        |             d           |
                #  way1 ...            ...        ... way2
                #        |               e         |
                #        |                \        |
                #        f-----------------g-------h
                a = external1[-2]
                b = common[0]
                c = external2[-2]
                d = common[1]
                e = common[-2]
                f = external1[1]
                g = common[-1]
                h = external2[1]

                data = [ SegmentedDetector.angle_abc(a,b,c),
                         SegmentedDetector.angle_abc(f,g,h),
                         SegmentedDetector.angle_abc(a,b,d),
                         SegmentedDetector.angle_abc(e,g,f),
                         SegmentedDetector.angle_abc(c,b,d),
                         SegmentedDetector.angle_abc(e,g,h)]

                data = [angle * 180 / math.pi for angle in data]
                data.extend([SegmentedDetector.diff_to_90(angle) for angle in data])

                # Compare common length ratio
                common_length = LineString(common).length
                external1_length = LineString(external1).length
                external2_length = LineString(external2).length
                ratio1 = common_length / external1_length
                ratio2 = common_length / external2_length
                data.extend([ratio1 + ratio2 / 2, min(ratio1, ratio2), max(ratio1, ratio2)])

                # Extended common part as they are with the cut on each side:
                common1_extd = [a] + common + [f]
                common2_extd = [c] + common + [h]
                # Consider extended ways, as they would be without the cut:
                external1_extd = [h] + external1 + [c]
                external2_extd = [f] + external2 + [a]

                external1_extd_angles, external2_extd_angles, common1_extd_angles, common2_extd_angles = \
                    [ numpy.array([SegmentedDetector.angle_abc(nodes[i-1], nodes[i], nodes[i+1]) * 180 / math.pi for i in xrange(1, len(nodes)-1)])
                      for nodes in external1_extd, external2_extd, common1_extd, common2_extd]

                data.extend(
                    [external1_extd_angles.mean(), external1_extd_angles.std(), external1_extd_angles.min(), external1_extd_angles.max(),
                     external2_extd_angles.mean(), external2_extd_angles.std(), external2_extd_angles.min(), external2_extd_angles.max(),
                     common1_extd_angles.mean() - external1_extd_angles.mean(),
                     common1_extd_angles.std(),
                     common1_extd_angles.min() - external1_extd_angles.min(),
                     common1_extd_angles.max() - external1_extd_angles.max(),
                     common2_extd_angles.mean() - external2_extd_angles.mean(),
                     common2_extd_angles.std(),
                     common2_extd_angles.min() - external2_extd_angles.min(),
                     common2_extd_angles.max() - external2_extd_angles.max()])

                external1_extd_angles, external2_extd_angles, common1_extd_angles, common2_extd_angles = \
                    [numpy.array([SegmentedDetector.diff_to_90(angle) for angle in angles]) for angles in
                        external1_extd_angles, external2_extd_angles, common1_extd_angles, common2_extd_angles ]

                data.extend(
                    [external1_extd_angles.mean(), external1_extd_angles.std(), external1_extd_angles.min(), external1_extd_angles.max(),
                     external2_extd_angles.mean(), external2_extd_angles.std(), external2_extd_angles.min(), external2_extd_angles.max(),
                     common1_extd_angles.mean() - external1_extd_angles.mean(),
                     common1_extd_angles.std(),
                     common1_extd_angles.min() - external1_extd_angles.min(),
                     common1_extd_angles.max() - external1_extd_angles.max(),
                     common2_extd_angles.mean() - external2_extd_angles.mean(),
                     common2_extd_angles.std(),
                     common2_extd_angles.min() - external2_extd_angles.min(),
                     common2_extd_angles.max() - external2_extd_angles.max()])

                result = data
        return result

    @staticmethod
    def diff_to_90(a):
        return abs(45-abs(45-(a%90)))


    @staticmethod
    def angle_abc(a,b,c):
        v1 = numpy.array([a[0]-b[0], a[1]-b[1]])
        v2 = numpy.array([c[0]-b[0], c[1]-b[1]])
        d = numpy.linalg.norm(v1) * numpy.linalg.norm(v2)
        if d == 0:
            return 0
        else:
            return numpy.arccos(numpy.clip(numpy.dot(v1, v2) / d, -1.0, 1.0))

    @staticmethod
    def get_external1_common_external2_ways(way1, way2):
        "return the part of way1 not common with way2, the common part, and the part of way2 not common with way1"
        assert(way1[-1] == way1[0]) # closed way
        assert(way2[-1] == way2[0]) # closed way
        way1 = way1[:-1]
        way2 = way2[:-1]
        previous_i = len(way1)-1
        for i in xrange(len(way1)):
            if way1[previous_i] not in way2 and way1[i] in way2:
                j = way2.index(way1[i])
                if (way2[(j + 1) % len(way2)] == way1[previous_i]) or \
                   (way2[(j - 1 + len(way2)) % len(way2)] == way1[(i+1) % len(way1)]):
                    # way2 is in reverse order
                    way2.reverse()
                    j = way2.index(way1[i])
                way1 = way1[i:] + way1[:i]
                way2 = way2[j:] + way2[:j]
                break
            previous_i = i
        i = 0
        while i<min(len(way1),len(way2)) and (way1[i] == way2[i]):
            i = i + 1
        if i==0:
           return way1+way1[0:1], [], way2 + way2[0:1]
        else:
           return way1[i-1:]+way1[0:1], way1[:i], way2[i-1:] + way2[0:1]


class Timer():
    def __init__(self, msg):
        self.start = timeit.default_timer()
        self.msg = msg
        print msg
    def __call__(self):
        return timeit.default_timer() - self.start
    def prnt(self):
        print self.msg + " => " + str(round(self(), 4)) +  " s"

