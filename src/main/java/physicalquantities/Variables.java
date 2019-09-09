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

public class Variables {

	private static Variables uniqueInstance;
	
	public static Variables getInstance() {
		/*if (uniqueInstance == null) {
			uniqueInstance = new Variables(waterSuction, temperature);
		}*/
		return uniqueInstance;
	}
	
	public static Variables getInstance(Map<Integer, Double> waterSuction) {
		if (uniqueInstance == null) {
			uniqueInstance = new Variables(waterSuction);
		}
		return uniqueInstance;
	}
	
	
	public static Map<Integer, Double> waterSuctions;
	public static Map<Integer, Double> thetasOld;
	public static Map<Integer, Double> thetas;
	public static Map<Integer, Double> thetasNew;
	public static Map<Integer, Double> saturationDegree;
	public static Map<Integer, Double> dThetas;
	public static Map<Integer, Double> dThetas1;
	public static Map<Integer, Double> dThetas2;
	public static Map<Integer, Double> thetas1;
	public static Map<Integer, Double> thetas2;
	public static Map<Integer, Double> kappas;
	public static Map<Integer, Double> volumes;
	public static Map<Integer, Double> volumesNew;
	public static Map<Integer, Double> darcyVelocities;
	public static Map<Integer, Double> darcyVelocitiesX;
	public static Map<Integer, Double> darcyVelocitiesZ;
	public static Map<Integer, Double> darcyVelocitiesCapillary;
	public static Map<Integer, Double> darcyVelocitiesGravity;
	public static Map<Integer, Double> poreVelocities;
	public static Map<Integer, Double> celerities;        // Rasmussen et al. 2000
	public static Map<Integer, Double> kinematicRatio;  // Rasmussen et al. 2000
	
	public static double errorVolume;
	public static double timeDelta;

	private Variables(Map<Integer, Double> waterSuction) {
		
		Variables.waterSuctions = new HashMap<Integer, Double>(waterSuction);
		Variables.thetasOld = new HashMap<Integer, Double>();
		Variables.thetas = new HashMap<Integer, Double>();
		Variables.thetasNew = new HashMap<Integer, Double>();
		Variables.saturationDegree = new HashMap<Integer, Double>();
		Variables.dThetas = new HashMap<Integer, Double>();
		Variables.dThetas1 = new HashMap<Integer, Double>();
		Variables.dThetas2 = new HashMap<Integer, Double>();
		Variables.thetas1 = new HashMap<Integer, Double>();
		Variables.thetas2 = new HashMap<Integer, Double>();
		Variables.kappas = new HashMap<Integer, Double>();
		Variables.volumes = new HashMap<Integer, Double>();
		Variables.volumesNew = new HashMap<Integer, Double>();
		Variables.darcyVelocities = new HashMap<Integer, Double>();
		Variables.darcyVelocitiesX = new HashMap<Integer, Double>();
		Variables.darcyVelocitiesZ = new HashMap<Integer, Double>();
		Variables.darcyVelocitiesCapillary = new HashMap<Integer, Double>();
		Variables.darcyVelocitiesGravity = new HashMap<Integer, Double>();
		Variables.poreVelocities = new HashMap<Integer, Double>();
		Variables.celerities = new HashMap<Integer, Double>();
		Variables.kinematicRatio = new HashMap<Integer, Double>();

	}

}
