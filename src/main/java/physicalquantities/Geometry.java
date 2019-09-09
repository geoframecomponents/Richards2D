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

package physicalquantities;

import java.util.HashMap;
import java.util.Map;

public class Geometry {

private static Geometry uniqueInstance;
	
	public static Geometry getInstance() {
		/*if (uniqueInstance == null) {
			uniqueInstance = new Geometry(elementsArea, edgesLenght, delta_j);
		}*/
		return uniqueInstance;
	}
	
	public static Geometry getInstance(Map<Integer, Double> elementsArea, Map<Integer, Double> edgesLenght, Map<Integer, Double> delta_j,
			Map<Integer, Double[]> edgeNormalVector, Map<Integer, Double[]> elementsCentroidsCoordinates, Map<Integer, Double[]> edgesCentroidsCoordinates) {
		if (uniqueInstance == null) {
			uniqueInstance = new Geometry(elementsArea, edgesLenght, delta_j, edgeNormalVector, elementsCentroidsCoordinates, edgesCentroidsCoordinates);
		}
		return uniqueInstance;
	}
	
	public static Map<Integer, Double> elementsArea;
	public static Map<Integer, Double> edgesLenght;
	public static Map<Integer, Double> delta_j;
	public static Map<Integer, Double[]> edgeNormalVector;
	public static Map<Integer, Double[]> elementsCentroidsCoordinates;
	public static Map<Integer, Double[]> edgesCentroidsCoordinates;


	
	private Geometry(Map<Integer, Double> elementsArea, Map<Integer, Double> edgesLenght, Map<Integer, Double> delta_j,
			 Map<Integer, Double[]> edgeNormalVector, Map<Integer, Double[]> elementsCentroidsCoordinates, Map<Integer, Double[]> edgesCentroidsCoordinates) {
		
		Geometry.elementsArea = new HashMap<Integer, Double>();
		Geometry.edgesLenght = new HashMap<Integer, Double>();
		Geometry.delta_j = new HashMap<Integer, Double>();
		Geometry.edgeNormalVector = new HashMap<Integer, Double[]>();
		Geometry.elementsCentroidsCoordinates = new HashMap<Integer, Double[]>();
		
		Geometry.elementsArea = elementsArea;
		Geometry.edgesLenght = edgesLenght;
		Geometry.delta_j = delta_j;
		Geometry.edgeNormalVector = edgeNormalVector;
		Geometry.elementsCentroidsCoordinates = elementsCentroidsCoordinates;
		Geometry.edgesCentroidsCoordinates = edgesCentroidsCoordinates;
		
	}

}
