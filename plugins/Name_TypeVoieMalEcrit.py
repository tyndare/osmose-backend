#-*- coding: utf-8 -*-

###########################################################################
##                                                                       ##
## Copyrights Etienne Chové <chove@crans.org> 2009                       ##
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

class plugin:
    
    only_for = ["fr"]
    
    err_702    = 5020
    err_702_fr = u"Type de voie mal écrit"
    err_702_en = u"Badly written way type"
    
    def init(self, logger):
        import re
        self.ReTests = {}        
        self.ReTests[( 0, u"Allée ")]     = re.compile(u"^([Aa][Ll][Ll]?[EÉée][Ee]?|[Aa][Ll][Ll]\.) .*$")
        self.ReTests[( 1, u"Boulevard ")] = re.compile(u"^[Bb]([Oo][Uu][Ll][Ll]?[Ee]?)?[Vv]?([Aa][Rr])?[Dd]\.? .*$")
        self.ReTests[( 2, u"Avenue ")]    = re.compile(u"^[Aa][Vv][Ee][Nn][Uu][Ee] .*$")
        #self.ReTests[( 3, u"École ")]     = re.compile(u"^[EÉée][Cc][Oo][Ll][Ee] .*$")
        self.ReTests[( 4, u"Chemin ")]    = re.compile(u"^[Cc][Hh][Ee][Mm][Ii][Nn] .*$")
        self.ReTests[( 5, u"Route ")]     = re.compile(u"^[Rr][Oo][Uu][Tt][Ee] .*$")
        self.ReTests[( 6, u"Esplanade ")] = re.compile(u"^[EÉée][Ss][Pp][Ll][Aa][Nn][Aa][Dd][Ee] .*$")
        self.ReTests[( 7, u"Rue ")]       = re.compile(u"^[Rr][Uu][Ee] .*$")
        self.ReTests = self.ReTests.items()
        
    def node(self, data, tags):
        if u"name" not in tags:
            return
        name = tags["name"]
        for test in self.ReTests:
            if test[1].match(name) and not name.startswith(test[0][1]):
                return [(702, test[0][0], {"fr": test[0][1].strip()})]
            
    def way(self, data, tags, nds):
        return self.node(data, tags)

    def relation(self, data, tags):
        return self.node(data, tags)
