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
 * A simple design factory to create an UnsaturatedHydraulicConductivity object.
 * @author Niccolo` Tubini
 */

public class UnsaturatedHydraulicConductivityFactory {
	
	
	
	/**
	 * Create a new SoilParametrization object.
	 * @param type name of the unsaturated hydraulic conductivity model
	 * @param modelSWRC soil water retention curve
	 * @param soilParametrization class containing the soil parameters 
	 * @return soilPar
	 */
	public UnsaturatedHydraulicConductivity createUnsaturatedHydraulicConductivity (String type, SoilWaterRetentionCurve modelSWRC) {

		UnsaturatedHydraulicConductivity unsatHydraulicModel = null;
		if(type.equalsIgnoreCase("Mualem Van Genuchten") || type.equalsIgnoreCase("MualemVanGenuchten")){
			unsatHydraulicModel = new MualemVanGenuchten(modelSWRC);
		}
//		else if (type.equalsIgnoreCase("Mualem Brooks Corey") || type.equalsIgnoreCase("MualemBrooksCorey")){
//			unsatHydraulicModel = new MualemBrooksCorey(modelSWRC);
//		}
		//else if(type.equalsIgnoreCase("Kosugi")){
		//	unsatHydraulicPar = new KosugiUnimodal(soilPar);
		//}
		//else if(type.equalsIgnoreCase("Romano")){
		//	unsatHydraulicPar = new Romano(soilPar);
		//}

		return unsatHydraulicModel;
		
		}

}
