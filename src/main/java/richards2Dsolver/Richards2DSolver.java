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
	double sumBoundaryFlow;
	int picardIteration;

	
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

	public Richards2DSolver(String soilHydraulicModel, String typeUHCModel, String matopType, int MAXITER_NEWT, double newtonTolerance, int picardIteration, double cgTolerance, boolean checkData) {
		
		this.picardIteration = picardIteration;
		soilWaterRetentionCurveFactory = new SoilWaterRetentionCurveFactory();
		soilWaterRetentionCurve = soilWaterRetentionCurveFactory.createSoilParametrization(soilHydraulicModel);

		unsaturatedHydraulicConductivityFactory = new UnsaturatedHydraulicConductivityFactory();
		unsaturatedHydraulicConductivity = unsaturatedHydraulicConductivityFactory.createUnsaturatedHydraulicConductivity(typeUHCModel, soilWaterRetentionCurve);

		matopFactory = new MatopFactory();
		matop = matopFactory.createMatop(matopType);

		nestedNewton = new NestedNewton(1, newtonTolerance, MAXITER_NEWT, soilWaterRetentionCurve, matop, cgTolerance);

		rhss = new HashMap<Integer,Double>();

		this.checkData = checkData;
	}


	public void solve(HashMap<Integer, double[]> inBC) {

		//		System.out.println("RICHARDS 2D] ");


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

		sumBoundaryFlow = 0.0; 

		/*
		 * Compute the right-hand side. The flux matrix entries are computed
		 * with the matop
		 * Note that the rhs contains the flux due to the gravitational potential
		 * 
		 * FIXME:
		 * - move the boundary conditions within the loop over edges for r == 0;
		 * - compute once for ever the gradient of z for each edge
		 */
		for(Integer element : Topology.s_i.keySet()) {
			rhss.put(element, Variables.volumes.get(element) );
		}
		
		for(int picard=0; picard<picardIteration; picard++) {
			System.out.println("\tPicard iteration: " +picard);
			/*
			 * Compute kappa at time level n
			 */
			for(Integer element : Topology.s_i.keySet()) {
				//Variables.kappas.put(element, Math.max(unsaturatedHydraulicConductivity.hydraulicConductivity(element), 0.000000000001) );
				//Variables.kappas.put(element, unsaturatedHydraulicConductivity.hydraulicConductivity(element) );
				Variables.kappas.put(element, Math.max(unsaturatedHydraulicConductivity.hydraulicConductivity(element), Math.ulp(1.0)) );


			}
//System.out.print(Math.ulp(1.0) +"  "+ Math.ulp(1E-15) +"  "+ Math.ulp(1E-21));
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
			 * 
			 * FIXME:
			 * - move the boundary conditions within the loop over edges for r == 0;
			 * - compute once for ever the gradient of z for each edge
			 */
			for(Integer element : Topology.s_i.keySet()) {
				rhss.put(element, Variables.volumes.get(element) );
			}


			for(Integer edge : Topology.l.keySet()) {
				/*
				 * Fluxes through edges within the domain
				 */
				double sideFlux = 0.0;

				if(Topology.r.get(edge)==0) {
					sideFlux = 0.0;  // depends on the type of the boundary condition and its value
					rhss.put(Topology.l.get(edge), rhss.get(Topology.l.get(edge))+sideFlux);
				} else {
					sideFlux = Variables.timeDelta*Math.max( Variables.kappas.get(Topology.l.get(edge)), Variables.kappas.get(Topology.r.get(edge)) )*( Geometry.elementsCentroidsCoordinates.get(Topology.r.get(edge))[1] - Geometry.elementsCentroidsCoordinates.get(Topology.l.get(edge))[1] )/Geometry.delta_j.get(edge)*Geometry.edgesLenght.get(edge);
					rhss.put( Topology.l.get(edge), rhss.get(Topology.l.get(edge))+sideFlux );
					rhss.put( Topology.r.get(edge), rhss.get(Topology.r.get(edge))-sideFlux );

				}
			}

			sumBoundaryFlow = 0.0;
			for(Integer edge : Topology.edgesBoundaryBCType.keySet()) {
				/*
				 * FIXME: boundary conditions
				 * if 0 no flux
				 * if 1 neumann 
				 */
				double sideFlux = 0.0;

				if(Topology.edgesBoundaryBCType.get(edge)==0) {
					sideFlux = 0.0;  // no flux
					rhss.put(Topology.l.get(edge), rhss.get(Topology.l.get(edge))-sideFlux);
				} else if(Topology.edgesBoundaryBCType.get(edge)==1) { // neumann
//					System.out.println(inBC.get(Topology.edgesBoundaryBCValue.get(edge))[0]);

					sideFlux = Variables.timeDelta*Geometry.edgesLenght.get(edge)*inBC.get(Topology.edgesBoundaryBCValue.get(edge))[0];  //0.00000667;

					rhss.put( Topology.l.get(edge), rhss.get(Topology.l.get(edge))+sideFlux );
					sumBoundaryFlow += sideFlux;

				} else if(Topology.edgesBoundaryBCType.get(edge)==2) { // dirichlet
//					System.out.println(inBC.get(Topology.edgesBoundaryBCValue.get(edge))[0]);
					double kappa =  Variables.kappas.get(Topology.l.get(edge)) * Geometry.edgesLenght.get(edge); //, unsaturatedHydraulicConductivity.hydraulicConductivity(0.0,Topology.l.get(edge)) 
					//				sideFlux = Variables.timeDelta * kappa * ( inBC.get(Topology.edgesBoundaryBCValue.get(edge))[0]/Geometry.delta_j.get(edge) + Geometry.edgeNormalVector.get(edge)[1] );// min(0.0,Geometry.edgeNormalVector.get(edge)[1])???
					sideFlux = Variables.timeDelta * kappa * ( inBC.get(Topology.edgesBoundaryBCValue.get(edge))[0]/Geometry.delta_j.get(edge) + ( Geometry.edgesCentroidsCoordinates.get(edge)[1] - Geometry.elementsCentroidsCoordinates.get(Topology.l.get(edge))[1] )/Geometry.delta_j.get(edge) );

					//				System.out.println("\tedge: " + edge);
					//				System.out.println("\t\tGeometry.edgesCentroidsCoordinates.get(edge)[1]: " + Geometry.edgesCentroidsCoordinates.get(edge)[1] + "\tGeometry.elementsCentroidsCoordinates.get(Topology.l.get(edge))[1]: " + Geometry.elementsCentroidsCoordinates.get(Topology.l.get(edge))[1] + "\tGeometry.delta_j.get(edge): " + Geometry.delta_j.get(edge));
					//				System.out.println("\t\tkappa: " + kappa + "\t :" + kappa * ( ( Geometry.edgesCentroidsCoordinates.get(edge)[1] - Geometry.elementsCentroidsCoordinates.get(Topology.l.get(edge))[1] )/Geometry.delta_j.get(edge) ));

					rhss.put( Topology.l.get(edge), rhss.get(Topology.l.get(edge))+sideFlux );

				} else if(Topology.edgesBoundaryBCType.get(edge)==3) { // free drainage
					//				System.out.println("\tfree drainage");
					double kappa =  Variables.kappas.get(Topology.l.get(edge)) * Geometry.edgesLenght.get(edge);
					/*
					 * FIXME: Math.min is necessary since the z component of normal vector can be positive and in that case
					 * the free drainage is 0
					 */
					sideFlux = Variables.timeDelta * kappa * (  Math.min(0.0,Geometry.edgeNormalVector.get(edge)[1]) );  // min(0.0,Geometry.edgeNormalVector.get(edge)[1])???
					rhss.put( Topology.l.get(edge), rhss.get(Topology.l.get(edge))+sideFlux );
					sumBoundaryFlow += sideFlux;
				}

				//sumBoundaryFlow += sideFlux;
			}


			if(checkData == true) { //change to true
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
			rhss.clear();
		} // close Picard iteration
		
		/*
		 * Update variables and compute the error
		 */
		/*
		 * Compute theta and saturation degree at time level n+1
		 */
		for(Integer element : Topology.s_i.keySet()) {
			Variables.thetas.put(element, soilWaterRetentionCurve.waterContent(element));
			Variables.saturationDegree.put(element, (Variables.thetas.get(element)-SoilParameters.thetaR[SoilParameters.elementsLabel.get(element)])/(SoilParameters.thetaS[SoilParameters.elementsLabel.get(element)]-SoilParameters.thetaR[SoilParameters.elementsLabel.get(element)]) );
		}

		if(checkData == true) {
			System.out.println("Theta at time level n+1:");
			for(Integer element : Variables.thetas.keySet()) {
				System.out.println("\t"+ element + "\t" + Variables.thetas.get(element));
			}
		}

		/*
		 * Compute water volumes at time level n+1
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

		/*
		 * Compute darcy fluxes at time level n+1
		 */
		double tmp;
		for(Integer edge : Topology.l.keySet()) {
			if(Topology.r.get(edge)==0) {
				Variables.darcyVelocities.put(edge,0.0);
				Variables.darcyVelocitiesX.put(edge,0.0*Geometry.edgeNormalVector.get(edge)[0]);
				Variables.darcyVelocitiesZ.put(edge,0.0*Geometry.edgeNormalVector.get(edge)[1]);
			} else {				
				double kappa = Math.max( Variables.kappas.get(Topology.l.get(edge)), Variables.kappas.get(Topology.r.get(edge)) )*Geometry.edgesLenght.get(edge);
				tmp = kappa * ( Variables.waterSuctions.get(Topology.r.get(edge)) + Geometry.elementsCentroidsCoordinates.get(Topology.r.get(edge))[1] -Variables.waterSuctions.get(Topology.l.get(edge)) - Geometry.elementsCentroidsCoordinates.get(Topology.l.get(edge))[1] )/Geometry.delta_j.get(edge);
				Variables.darcyVelocities.put(edge,tmp);
				//				Variables.darcyVelocitiesX.put(edge,tmp*Math.abs(Geometry.edgeNormalVector.get(edge)[0]));
				//				Variables.darcyVelocitiesZ.put(edge,tmp*Math.abs(Geometry.edgeNormalVector.get(edge)[1]));
				if(Math.abs(Geometry.elementsCentroidsCoordinates.get(Topology.r.get(edge))[0]-Geometry.elementsCentroidsCoordinates.get(Topology.l.get(edge))[0]) < 0.0000001) {
					Variables.darcyVelocitiesX.put(edge,0.0);
				} else {
					Variables.darcyVelocitiesX.put(edge, -kappa*Math.abs(Geometry.edgeNormalVector.get(edge)[0])*( Variables.waterSuctions.get(Topology.r.get(edge)) + Geometry.elementsCentroidsCoordinates.get(Topology.r.get(edge))[1] -Variables.waterSuctions.get(Topology.l.get(edge)) - Geometry.elementsCentroidsCoordinates.get(Topology.l.get(edge))[1] )/(Geometry.elementsCentroidsCoordinates.get(Topology.r.get(edge))[0]-Geometry.elementsCentroidsCoordinates.get(Topology.l.get(edge))[0]) );
				}

				if(Math.abs(Geometry.elementsCentroidsCoordinates.get(Topology.r.get(edge))[1]-Geometry.elementsCentroidsCoordinates.get(Topology.l.get(edge))[1]) < 0.0000001) {
					Variables.darcyVelocitiesZ.put(edge,0.0);
				} else {
					Variables.darcyVelocitiesZ.put(edge, -kappa*Math.abs(Geometry.edgeNormalVector.get(edge)[1])*( Variables.waterSuctions.get(Topology.r.get(edge)) + Geometry.elementsCentroidsCoordinates.get(Topology.r.get(edge))[1] -Variables.waterSuctions.get(Topology.l.get(edge)) - Geometry.elementsCentroidsCoordinates.get(Topology.l.get(edge))[1] )/(Geometry.elementsCentroidsCoordinates.get(Topology.r.get(edge))[1]-Geometry.elementsCentroidsCoordinates.get(Topology.l.get(edge))[1]) );
				}
			}
		}


		/*
		 * Update border fluxes with dirichlet 
		 */
		for(Integer edge : Topology.edgesBoundaryBCType.keySet()) {
			if(Topology.edgesBoundaryBCType.get(edge)==2) { // 
				double kappa =  Variables.kappas.get(Topology.l.get(edge)) * Geometry.edgesLenght.get(edge); //, unsaturatedHydraulicConductivity.hydraulicConductivity(0.0,Topology.l.get(edge)) 
				//sideFlux = Variables.timeDelta * kappa * ( 0.0 - Geometry.elementsCentroidsCoordinates.get(Topology.l.get(edge))[1] )/Geometry.delta_j.get(edge);
				sumBoundaryFlow += Variables.timeDelta * kappa * ( (inBC.get(Topology.edgesBoundaryBCValue.get(edge))[0]-Variables.waterSuctions.get(Topology.l.get(edge)))/Geometry.delta_j.get(edge) + Geometry.edgeNormalVector.get(edge)[1] );
			}
		}
		double volume = 0.0;
		double volumeNew = 0.0;
		for(Integer element : Topology.s_i.keySet()) {
			volume += Variables.volumes.get(element);
			volumeNew += Variables.volumesNew.get(element);
		}

		double errorVolume = volumeNew - volume - (sumBoundaryFlow + 0.0);
		System.out.println("\tERROR VOLUME : " + errorVolume);

	}

}
