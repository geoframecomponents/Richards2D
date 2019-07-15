/*
 * GNU GPL v3 License
 *
 * Copyright 2017  Niccolo` Tubini
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
package rehology;


/**
 * A simple design factory to create a SoilParametrization object.
 * @author Niccolo` Tubini
 */

public class SoilWaterRetentionCurveFactory {
	
	
	
	/**
	 * Create a new SoilParametrization object.
	 * @param type name of the SWRC model
	 * @return soilPar
	 */
	public SoilWaterRetentionCurve createSoilParametrization (String type) {

		SoilWaterRetentionCurve soilPar = null;
//		if(type.equalsIgnoreCase("BrooksCorey") || type.equalsIgnoreCase("Brooks Corey")){
//			soilPar = new BrooksCorey();
//		}
		if(type.equalsIgnoreCase("VanGenuchten") || type.equalsIgnoreCase("Van Genuchten")){
			soilPar = new VanGenuchten();
		}
//		else if(type.equalsIgnoreCase("Kosugi")){
//			soilPar = new KosugiUnimodal();
//		}
//		else if(type.equalsIgnoreCase("Romano")){
//			soilPar = new Romano();
//		}
//		else if(type.equalsIgnoreCase("VanGenuchtenB") || type.equalsIgnoreCase("Van Genuchten B")){
//			soilPar = new VanGenuchtenBachmann();
//		}
		return soilPar;
		}

}
