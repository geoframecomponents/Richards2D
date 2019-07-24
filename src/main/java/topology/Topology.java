/*
 * GNU GPL v3 License
 *
 * Copyright 2019 Niccolo` Tubini
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package topology;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

public class Topology {

private static Topology uniqueInstance;
	
	public static Topology getInstance() {
		/*if (uniqueInstance == null) {
			uniqueInstance = new Geometry(elementsArea, edgesLenght, delta_j);
		}*/
		return uniqueInstance;
	}
	
	public static Topology getInstance(Map<Integer, Integer> l, Map<Integer, Integer> r, Map<Integer, Integer> edgesBoundaryBCType,
			 							Map<Integer, Integer> edgesBoundaryBCValue, Map<Integer, ArrayList<Integer>> s_i) {
		if (uniqueInstance == null) {
			uniqueInstance = new Topology(l, r, edgesBoundaryBCType, edgesBoundaryBCValue, s_i);
		}
		return uniqueInstance;
	}
	
	public static Map<Integer, Integer> l;
	public static Map<Integer, Integer> r;
	public static Map<Integer, Integer> edgesBoundaryBCType;
	public static Map<Integer, Integer> edgesBoundaryBCValue;
	public static Map<Integer, ArrayList<Integer>> s_i;

	
	private Topology(Map<Integer, Integer> l, Map<Integer, Integer> r, Map<Integer, Integer> edgesBoundaryBCType,
			 			Map<Integer, Integer> edgesBoundaryBCValue, Map<Integer, ArrayList<Integer>> s_i) {
		
		Topology.l = new HashMap<Integer, Integer>(l);
		Topology.r = new HashMap<Integer, Integer>(r);
		Topology.edgesBoundaryBCType = new HashMap<Integer, Integer>(edgesBoundaryBCType);
		Topology.edgesBoundaryBCValue = new HashMap<Integer, Integer>(edgesBoundaryBCValue);
		Topology.s_i = new HashMap<Integer, ArrayList<Integer>>(s_i);
				
	}

}
