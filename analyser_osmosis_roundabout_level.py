#!/usr/bin/env python
#-*- coding: utf-8 -*-

###########################################################################
##                                                                       ##
## Copyrights Frédéric Rodrigo 2011                                      ##
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

sql10 = """
CREATE OR REPLACE FUNCTION level(highway varchar) RETURNS int AS $$
DECLARE BEGIN
    RETURN CASE
        WHEN highway = 'motorway' THEN 7
        WHEN highway = 'trunk' THEN 7
        WHEN highway = 'primary' THEN 6
        WHEN highway = 'secondary' THEN 5
        WHEN highway = 'tertiary' THEN 4
        WHEN highway = 'unclassified' THEN 3
        WHEN highway = 'residential' THEN 3
        WHEN highway = 'service' THEN 2
        WHEN highway = 'road' THEN 1
        ELSE 0
    END;
END
$$ LANGUAGE plpgsql;

DROP VIEW IF EXISTS roundabout CASCADE;
CREATE VIEW roundabout AS
SELECT
    id,
    level(tags->'highway') AS level,
    tags->'highway' AS highway,
    linestring
FROM
    ways
WHERE
    tags?'junction' AND
    tags->'junction' = 'roundabout' AND
    tags?'highway' AND
    nodes[1] = nodes[array_length(nodes,1)]
;
"""

sql11 = """
SELECT
    roundabout.id,
    ST_AsText(way_locate(roundabout.linestring)),
    roundabout.level
FROM
    roundabout
    JOIN way_nodes AS wn1 ON
        roundabout.id = wn1.way_id
    JOIN way_nodes AS wn2 ON
        wn1.node_id = wn2.node_id AND
        roundabout.id != wn2.way_id
    JOIN ways ON
        wn2.way_id = ways.id
WHERE
    ways.tags?'highway'
GROUP BY
    roundabout.id,
    roundabout.level,
    roundabout.highway,
    roundabout.linestring
HAVING
    MAX(level(tags->'highway')) < 7 AND -- doesn't force motorway or trunk roundabout as local trafic may pass through
    MAX(level(tags->'highway')) != roundabout.level
;
"""

sql20 = """
CREATE OR REPLACE FUNCTION other_end(n_id bigint, nodes bigint[]) RETURNS bigint AS $$
DECLARE BEGIN
    IF nodes[1] = n_id THEN
        RETURN nodes[array_length(nodes,1)];
    END IF;
    IF nodes[array_length(nodes,1)] = n_id THEN
        RETURN nodes[1];
    END IF;
    RETURN NULL;
END
$$ LANGUAGE plpgsql;

DROP TABLE IF EXISTS roundabout_acces;
CREATE TEMP TABLE roundabout_acces AS
SELECT
    roundabout.id AS ra_id,
    ways.id AS a_id,
    other_end(wn1.node_id, ways.nodes) AS n_id,
    ways.linestring,
    (ways.tags?'oneway' AND ways.tags->'oneway' IN ('yes', 'true', '-1', '1')) AS oneway
FROM
    roundabout
    JOIN way_nodes AS wn1 ON
        roundabout.id = wn1.way_id
    JOIN way_nodes AS wn2 ON
        roundabout.id != wn2.way_id AND
        wn1.node_id = wn2.node_id
    JOIN ways ON
        wn2.way_id = ways.id
WHERE
    ways.tags?'highway' AND
    ways.tags->'highway' IN ('primary', 'secondary', 'tertiary', 'unclassified', 'residential', 'road') AND
    array_length(ways.nodes,1) <= 4
;

CREATE INDEX roundabout_acces_idx ON roundabout_acces(ra_id);
"""

sql21 = """
SELECT
    ra1.a_id,
    ST_AsText(way_locate(ra1.linestring))
FROM
    roundabout_acces AS ra1,
    roundabout_acces AS ra2
WHERE
    ra1.ra_id = ra2.ra_id AND
    ra1.a_id != ra2.a_id AND
    ra1.n_id = ra2.n_id AND
    NOT ra1.oneway
GROUP BY
    ra1.a_id,
    ra1.linestring
;
"""

sql30 = """
SELECT
    junction.id AS junction_id,
    ST_AsText(ST_Centroid(ST_Intersection(w1.linestring, w2.linestring))) -- Centroid beacause can be any geom
FROM
    ways AS junction
    JOIN ways AS w1 ON
        junction.linestring && w1.linestring AND
        w1.tags?'highway' AND
        w1.id != junction.id
    JOIN ways AS w2 ON
        junction.linestring && w2.linestring AND
        w2.tags?'highway' AND
        w2.id != junction.id
WHERE
    array_length(junction.nodes, 1) > 3 AND
    junction.nodes[1] = junction.nodes[array_length(junction.nodes, 1)] AND
    w1.id != w2.id AND
    junction.tags?'junction' AND
    junction.tags->'junction' = 'roundabout' AND
    w1.linestring && w2.linestring AND
    (select array_agg(e) from (SELECT unnest(junction.nodes) INTERSECT SELECT ends(w1.nodes)) AS dt(e)) =
    (select array_agg(e) from (SELECT unnest(junction.nodes) INTERSECT SELECT ends(w2.nodes)) AS dt(e))
;
"""


class Analyser_Osmosis_Roundabout_Level(Analyser_Osmosis):

    def __init__(self, config, logger = None):
        Analyser_Osmosis.__init__(self, config, logger)
        self.classs[1] = {"item":"3010", "level": 2, "tag": ["highway", "roundabout"], "desc":{"fr":"Mauvais highway sur roundabout", "en":"Wrong highway on roundabout"} }
        self.classs[2] = {"item":"2030", "level": 2, "tag": ["highway", "roundabout"], "desc":{"fr":"oneway manquant sur insertion Rond-Point", "en":"Missing oneway"} }
        self.classs[3] = {"item":"3010", "level": 2, "tag": ["highway", "roundabout"], "desc":{"fr":"Raccourci sur rond-point", "en":"Roundabout shortcut"} }

    def analyser_osmosis(self):
        self.run(sql10)
        self.run(sql11, lambda res: {"class":1, "subclass":res[2], "data":[self.way_full, self.positionAsText]} )
        self.run(sql20)
        self.run(sql21, lambda res: {"class":2, "data":[self.way_full, self.positionAsText]} )
        self.run(sql30, lambda res: {"class":3, "data":[self.way_full, self.positionAsText]} )
