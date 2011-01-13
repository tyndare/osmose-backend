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
    
    err_105    = 1020
    err_105_fr = u"Way sans nœud"
    err_105_en = u"Way without node"
    
    err_106    = 1030
    err_106_fr = u"Way avec un seul nœud"
    err_106_en = u"Way with a single node"
    
    def way(self, data, tags, nds):
        if len(nds) >= 2:
            return
        way = self.father.WayGet(data[u"id"])
        if not way:
            return
        realnb = len(way[u"nd"])
        if realnb == 0:
            return [(105, 0, {})]
        if realnb == 1:
            return [(106, 0, {})]
