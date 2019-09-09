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
import java.util.Map;

import oms3.annotations.*;
import physicalquantities.*;
import topology.Topology;

public class CallRichards2DSolver {

	@In
	public Map<Integer, Integer> l;

	@In
	public Map<Integer, Integer> r;

	@In
	public Map<Integer, ArrayList<Integer>> s_i;

	@In
	public Map<Integer, Double> elementsArea;

	@In
	public Map<Integer, Double> edgesLenght;

	@In
	public Map<Integer, Double> delta_j;

	@In
	public Map<Integer, Double[]> edgeNormalVector;
	
	@In
	public Map<Integer, Double[]> elementsCentroidsCoordinates;
	
	@In
	public Map<Integer, Double[]> edgesCentroidsCoordinates;
	
	@In
	public Map<Integer, Integer> elementsLabel;

	@In
	public Map<Integer, Integer> edgesBoundaryBCType;
	
	@In
	public Map<Integer, Integer> edgesBoundaryBCValue;

	@In
	public double[] alphaSpecificStorage;

	@In
	public double[] betaSpecificStorage;

	@In
	public double[] ks;

	@In
	public double[] par1SWRC;

	@In
	public double[] par2SWRC;

	@In
	public double[] par3SWRC;

	@In
	public double[] par4SWRC;

	@In
	public double[] par5SWRC;

	@In
	public double[] psiStar1;

	@In
	public double[] psiStar2;

	@In
	public double[] psiStar3;

	@In
	public double[] thetaR;

	@In
	public double[] thetaS;

	//@In
	public Map<Integer, Double> psi;

	@In
	public String soilHydraulicModel;

	@In
	public String typeUHCModel;
	
	@In
	public String typeMatop;
	
	@In
	public String initialConditionType;

	@In
	public double newtonTolerance = 0.000000001;
	
	@In
	public int MAXITER_NEWT = 20;
	
	@In
	public int picardIteration =1;
	
	@In
	public double cgTolerance = 0.000000001;

	@In
	public boolean checkData = false;
	
	@Description("Time amount at every time-loop")
	@In
	@Unit ("s")
	public double tTimestep = 1.0;

	@Description("Time step of integration")
	@In
	@Unit ("s")
	public double timeDelta = 1.0;
	
	@Description("Position of the seepage")
	@In
	@Unit ("m")
	public double zSeePage = 0.0;
	
	// BOUNDARY CONDITIONS

	@Description("The HashMap with the time series of the boundary condition at the top of soil column")
	@In
	@Unit ("m")
	public HashMap<Integer, double[]> inBC;
	
	@Description("The current date of the simulation.")
	@In
	public String inCurrentDate;
	
	@Description("ArrayList of variable to be stored in the buffer writer")
	@Out
	public ArrayList<double[]> outputToBuffer;
	
	//////////////////////////////////////
	//////////////////////////////////////

	Variables variables;
	SoilParameters soilParameters;
	Geometry geometry;
	Topology topology;

	double[] tmpElement;
	double[] tmpEdge;
	
	int step = 0;

	///////////////////

	Richards2DSolver richardsSolver;

	@Execute
	public void solve() {

		System.out.println("RICHARDS 2D " + inCurrentDate);

		if(step==0){
			
			psi = new HashMap<Integer, Double>();
			/*
			 * Compute the initial condition
			 */
			if( initialConditionType.equalsIgnoreCase("hydrostatic") ) {
				for(Integer i : elementsCentroidsCoordinates.keySet()) {
					psi.put(i, zSeePage-elementsCentroidsCoordinates.get(i)[1]);
				}
			} else if (initialConditionType.equalsIgnoreCase("constant") ) {
				for(Integer i : elementsCentroidsCoordinates.keySet()) {
					psi.put(i,zSeePage);
				}
			} else {
				System.out.println("Error: initial condition not valid.");
			}

			variables = Variables.getInstance(psi);
			Variables.timeDelta = timeDelta;
			soilParameters = SoilParameters.getInstance(elementsLabel, alphaSpecificStorage, betaSpecificStorage, ks, 
					par1SWRC, par2SWRC, par3SWRC, par4SWRC, par5SWRC, psiStar1, psiStar2, psiStar3, thetaR, thetaS);
			geometry = Geometry.getInstance(elementsArea, edgesLenght, delta_j, edgeNormalVector, elementsCentroidsCoordinates, edgesCentroidsCoordinates);
			topology = Topology.getInstance(l, r, edgesBoundaryBCType, edgesBoundaryBCValue, s_i);

			richardsSolver = new Richards2DSolver(soilHydraulicModel, typeUHCModel, typeMatop, MAXITER_NEWT, newtonTolerance, picardIteration, cgTolerance, checkData);
			
			tmpElement = new double[Topology.s_i.size()+1];
			tmpEdge = new double[Topology.r.size()+1];
			
			outputToBuffer= new ArrayList<double[]>();
		}

		
		outputToBuffer.clear();
		
		double sumTimeDelta = 0;

		//double volume = 0.0;
		//double volumeNew = 0.0;
		while(sumTimeDelta < tTimestep) {

			if(sumTimeDelta + Variables.timeDelta>tTimestep) {
				Variables.timeDelta = tTimestep - sumTimeDelta;
			}
			sumTimeDelta = sumTimeDelta + Variables.timeDelta;
			
			richardsSolver.solve(inBC);
			
		}
		
		for(Integer i : psi.keySet()) {
			tmpElement[i] = psi.get(i);
		}
		outputToBuffer.add(tmpElement.clone());
		
		for(Integer i : Variables.waterSuctions.keySet()) {
			tmpElement[i] = Variables.waterSuctions.get(i);
		}
		outputToBuffer.add(tmpElement.clone());
		
		for(Integer i : Variables.thetas.keySet()) {
			tmpElement[i] = Variables.thetas.get(i);
		}
		outputToBuffer.add(tmpElement.clone());
		
		for(Integer i : Variables.saturationDegree.keySet()) {
			tmpElement[i] = Variables.saturationDegree.get(i);
		}
		outputToBuffer.add(tmpElement.clone());
		
		for(Integer i : Variables.darcyVelocities.keySet()) {
			tmpEdge[i] = Variables.darcyVelocities.get(i);
		}
		outputToBuffer.add(tmpEdge.clone());
		
		for(Integer i : Variables.darcyVelocitiesX.keySet()) {
			tmpEdge[i] = Variables.darcyVelocitiesX.get(i);
		}
		outputToBuffer.add(tmpEdge.clone());
		
		for(Integer i : Variables.darcyVelocitiesZ.keySet()) {
			tmpEdge[i] = Variables.darcyVelocitiesZ.get(i);
		}
		outputToBuffer.add(tmpEdge.clone());
				
		
		
		step ++;
		
		System.out.println("\n\n");
		
	}
}
