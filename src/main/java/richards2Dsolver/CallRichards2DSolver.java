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
	public Map<Integer, Integer> elementsLabel;

	@In
	public Map<Integer, Integer> edgesLabel;

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

	@In
	public Map<Integer, Double> psi;

	@In
	public String soilHydraulicModel;

	@In
	public String typeUHCModel;
	
	@In
	public String typeMatop;


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
	//////////////////////////////////////
	//////////////////////////////////////

	Variables variables;
	SoilParameters soilParameters;
	Geometry geometry;
	Topology topology;


	int step = 0;

	///////////////////

	Richards2DSolver richardsSolver;

	@Execute
	public void solve() {

		System.out.println("RICHARDS 2D ");

		if(step==0){
			variables = Variables.getInstance(psi);
			Variables.timeDelta = timeDelta;
			soilParameters = SoilParameters.getInstance(elementsLabel, alphaSpecificStorage, betaSpecificStorage, ks, 
					par1SWRC, par2SWRC, par3SWRC, par4SWRC, par5SWRC, psiStar1, psiStar2, psiStar3, thetaR, thetaS);
			geometry = Geometry.getInstance(elementsArea, edgesLenght, delta_j, edgeNormalVector, elementsCentroidsCoordinates);
			topology = Topology.getInstance(l, r, edgesLabel, s_i);

			richardsSolver = new Richards2DSolver(soilHydraulicModel, typeUHCModel, typeMatop, checkData);
		}

		double sumTimeDelta = 0;

		//double volume = 0.0;
		//double volumeNew = 0.0;
		while(sumTimeDelta < tTimestep) {

			if(sumTimeDelta + Variables.timeDelta>tTimestep) {
				Variables.timeDelta = tTimestep - sumTimeDelta;
			}
			sumTimeDelta = sumTimeDelta + Variables.timeDelta;
			
			richardsSolver.solve(0.0);
			
		}
		
		System.out.println("\n\nSOLUTION:");
		for(Integer element : Topology.s_i.keySet()) {
			System.out.println("\t" + element + "\t" + Variables.waterSuctions.get(element) );
		}
		System.out.println("\n\n");
	}
}
