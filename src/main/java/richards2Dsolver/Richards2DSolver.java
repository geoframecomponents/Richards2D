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

package richards2Dsolver;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import oms3.annotations.Description;
import physicalquantities.*;
import rehology.*;
import topology.Topology;

public class Richards2DSolver {

	Map<Integer, Double> rhss;

	double sigma;
	double tmp;
	double newtonTolerance = 0.00000000001;

	@Description("Object dealing with SWRC model")
	SoilWaterRetentionCurve soilWaterRetentionCurve;
	SoilWaterRetentionCurveFactory soilWaterRetentionCurveFactory;
	//SoilWaterRetentionCurveTemperatureFactory soilWaterRetentionCurveTemperatureFactory;

	@Description("Object dealing with unsaturated hydraulic conductivity model")
	UnsaturatedHydraulicConductivity unsaturatedHydraulicConductivity;
	UnsaturatedHydraulicConductivityFactory unsaturatedHydraulicConductivityFactory;
	//UnsaturatedHydraulicConductivityTemperatureFactory unsaturatedHydraulicConductivityTemperatureFactory;

	/*
	 * FIXME:
	 * mancano le classi per le boundary condition
	 */

	/*
	 * FIXME:
	 * classe per calcolare la conducibilita` idraulica tra volumi adiacenti? max min mean
	 */

	@Description("Object to perform the nested Newton algortithm")
	NestedNewton nestedNewton;

	@Description("Object to compute the product between the flux matrix and the variable vector")
	Matop matop;
	MatopFactory matopFactory;

	boolean checkData;

	public Richards2DSolver(String soilHydraulicModel, String typeUHCModel, String matopType, boolean checkData) {

		soilWaterRetentionCurveFactory = new SoilWaterRetentionCurveFactory();
		soilWaterRetentionCurve = soilWaterRetentionCurveFactory.createSoilParametrization(soilHydraulicModel);

		unsaturatedHydraulicConductivityFactory = new UnsaturatedHydraulicConductivityFactory();
		unsaturatedHydraulicConductivity = unsaturatedHydraulicConductivityFactory.createUnsaturatedHydraulicConductivity(typeUHCModel, soilWaterRetentionCurve);

		matopFactory = new MatopFactory();
		matop = matopFactory.createMatop(matopType);

		nestedNewton = new NestedNewton(1, newtonTolerance, 10, soilWaterRetentionCurve, matop);

		rhss = new HashMap<Integer,Double>();

		this.checkData = checkData;
	}


	public void solve(double timeDelta) {

		System.out.println("RICHARDS 2D] ");


		/*
		 * Compute theta at time level n
		 */
		for(Integer element : Topology.s_i.keySet()) {
			Variables.thetasOld.put(element, soilWaterRetentionCurve.waterContent(element));
		}

		if(checkData == true) {
			System.out.println("Theta at time level n:");
			for(Integer element : Variables.thetasOld.keySet()) {
				System.out.println("\t"+ element + "\t" + Variables.thetasOld.get(element));
			}
		}


		/*
		 * Compute water volumes at time level n
		 */
		for(Integer element : Topology.s_i.keySet()) {
			Variables.volumes.put(element, Variables.thetasOld.get(element)*Geometry.elementsArea.get(element));
		}

		if(checkData == true) {
			System.out.println("Volume at time level n:");
			for(Integer element : Variables.volumes.keySet()) {
				System.out.println("\t"+ element + "\t" + Variables.volumes.get(element) + "\t" + Geometry.elementsArea.get(element));
			}
		}

		
		for(int picard=0; picard<1; picard++) {

			/*
			 * Compute kappa at time level n
			 */
			for(Integer element : Topology.s_i.keySet()) {
				Variables.kappas.put(element, unsaturatedHydraulicConductivity.hydraulicConductivity(element));
			}

			if(checkData == true) {
				System.out.println("Kappa at time level n:");
				for(Integer element : Variables.kappas.keySet()) {
					System.out.println("\t"+ element + "\t" + Variables.kappas.get(element));
				}
			}


			/*
			 * Compute the right-hand side. The flux matrix entries are computed
			 * with the matop
			 * Note that the rhs contains the flux due to the gravitational potential
			 */
			for(Integer element : Topology.s_i.keySet()) {
				rhss.put(element, Variables.volumes.get(element) );
			}

			for(Integer edge : Topology.l.keySet()) {
				/*
				 * FIXME: everywhere NO FLUX BOUNDARY CONDITION
				 */
				double sideFlux = 0.0;

				if(Topology.r.get(edge)==0) {
					sideFlux = 0.0;  // depends on the type of the boundary condition and its value
					rhss.put(Topology.l.get(edge), rhss.get(Topology.l.get(edge))-sideFlux);
				} else {
					sideFlux = Variables.timeDelta*Math.max( Variables.kappas.get(Topology.l.get(edge)), Variables.kappas.get(Topology.r.get(edge)) )*( Geometry.elementsCentroidsCoordinates.get(Topology.r.get(edge))[1] - Geometry.elementsCentroidsCoordinates.get(Topology.l.get(edge))[1] )/Geometry.delta_j.get(edge)*Geometry.edgesLenght.get(edge);

					rhss.put( Topology.l.get(edge), rhss.get(Topology.l.get(edge))+sideFlux );
					rhss.put( Topology.r.get(edge), rhss.get(Topology.r.get(edge))-sideFlux );

				}
			}


			if(checkData == true) {
				System.out.println("rhss :");
				for(Integer element : Topology.s_i.keySet()) {
					System.out.println("\t"+ element + "\t" + rhss.get(element));
				}
			}

			/*
			 * Linearize and solve the system:
			 * - nested newton algorithm (Casulli and Zanolli, 2010)
			 * - conjugate gradient method ( Shewchuk, 1994)
			 */

			// Initial guess for the outer iteration
			nestedNewton.set(rhss);
			nestedNewton.solver();


			/*
			 * Update variables and compute the error
			 */
			/*
			 * Compute theta at time level n+1
			 */
			for(Integer element : Topology.s_i.keySet()) {
				Variables.thetas.put(element, soilWaterRetentionCurve.waterContent(element));
			}

			if(checkData == true) {
				System.out.println("Theta at time level n+1:");
				for(Integer element : Variables.thetas.keySet()) {
					System.out.println("\t"+ element + "\t" + Variables.thetas.get(element));
				}
			}

			/*
			 * Compute water volumes at time level n
			 */
			for(Integer element : Topology.s_i.keySet()) {
				Variables.volumesNew.put(element, Variables.thetas.get(element)*Geometry.elementsArea.get(element));
			}

			if(checkData == true) {
				System.out.println("Volume at time level n + 1:");
				for(Integer element : Variables.volumes.keySet()) {
					System.out.println("\t"+ element + "\t" + Variables.volumesNew.get(element) + "\t" + Geometry.elementsArea.get(element));
				}
			}

			double volume = 0.0;
			double volumeNew = 0.0;
			for(Integer element : Topology.s_i.keySet()) {
				volume += Variables.volumes.get(element);
				volumeNew += Variables.volumesNew.get(element);
			}

			double errorVolume = volumeNew - volume - timeDelta*(0.0 + 0.0);
			System.out.println("ERROR VOLUME : " + errorVolume);
			//			for(Integer element : Topology.s_i.keySet()) {
			//				System.out.println("\t" + element + "\t" + Variables.waterSuctions.get(element) );
			//			}
			//			System.out.println("\n\n");
		}
		//} // close picard iteration
	}

}
