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

import java.util.HashMap;
import java.util.Map;

import org.junit.Test;

import generatemesh.GenerateTriangularMesh;
import meshtopology.TopologyTriangularMesh;
import readtriangularization.Readmsh;

public class TestCallRichardsSolver {

	@Test
	public void Test() throws Exception {

		//String fileName = "resources/input/square_with_subdomain_100.msh";
		String fileName = "resources/input/square22_1.msh";

		String splitter = " ";

		Readmsh reader = new Readmsh();
		reader.fileName = fileName;
		reader.splitter = splitter;
		reader.checkData = false;
		reader.process();
		
		//GenerateMeshTriangles generateMesh = new GenerateMeshTriangles();
		GenerateTriangularMesh generateMesh = new GenerateTriangularMesh();
		generateMesh.verticesCoordinates = reader.verticesCoordinates;
		generateMesh.elementsVertices = reader.elementsVertices;
		generateMesh.borderEdgesVertices = reader.borderEdgesVertices;
		generateMesh.borderEdgesLabel = reader.borderEdgesLabel;
		generateMesh.checkData = false;
		generateMesh.geometryType = "EuclideanCartesian";
		generateMesh.process();
		
		double[] alphaSpecificStorage = new double[] {0.0, 0.0};
		double[] betaSpecificStorage = new double[] {0.0, 0.0};
		double[] ks = new double[] {0.28, 0.28};//{0.0000028, 0.0000123};
		double[] par1SWRC = new double[] {1.56,1.56};
		double[] par2SWRC = new double[] {3.6,3.6};
		double[] par3SWRC = null;
		double[] par4SWRC = null;
		double[] par5SWRC = null;
		double[] psiStar1 = new double[] {-0.144,-0.144};
		double[] psiStar2 = null;
		double[] psiStar3 = null;
		double[] thetaS = new double[] {0.43,0.43};
		double[] thetaR = new double[] {0.078,0.078};
		
		/*
		 * Compute the initial condition
		 */
		Map<Integer, Double> psi = new HashMap<Integer, Double>();
		for(Integer i : generateMesh.elementsCentroidsCoordinates.keySet()) {
			//psi.put(i, -generateMesh.elementsCentroidsCoordinates.get(i)[1]);
			psi.put(i,-0.5);
		}
		//psi.put(6, 0.0);
		//psi.put(5, 0.0);
		//psi.put(8, 0.0);
		//psi.put(7, 0.0);
		System.out.println("Initial condition:");
		for(Integer i : psi.keySet()) {
			System.out.println("\telement:" + i + " psi = " + psi.get(i));
		}
		
		CallRichards2DSolver solver = new CallRichards2DSolver();
		solver.l = generateMesh.l;
		solver.r = generateMesh.r;
		solver.s_i = generateMesh.s_i;		
		solver.elementsArea = generateMesh.elementsArea;
		solver.edgesLenght = generateMesh.edgesLength;
		solver.delta_j = generateMesh.delta_j;
		solver.edgeNormalVector = generateMesh.edgeNormalVector;
		solver.elementsCentroidsCoordinates = generateMesh.elementsCentroidsCoordinates;
		solver.alphaSpecificStorage = alphaSpecificStorage;
		solver.betaSpecificStorage = betaSpecificStorage;
		solver.ks = ks;
		solver.par1SWRC = par1SWRC;
		solver.par2SWRC = par2SWRC;
		solver.par3SWRC = par3SWRC;
		solver.par4SWRC = par4SWRC;
		solver.par5SWRC = par5SWRC;
		solver.psiStar1 = psiStar1;
		solver.psiStar2 = psiStar2;
		solver.psiStar3 = psiStar3;
		solver.thetaS = thetaS;
		solver.thetaR = thetaR;
		solver.elementsLabel = reader.elementsLabel;
		solver.edgesLabel = reader.borderEdgesLabel;
		solver.psi = psi;
		solver.soilHydraulicModel = "VanGenuchten";
		solver.typeUHCModel = "MualemVanGenuchten";
		solver.typeMatop = "2DRichards";
		solver.checkData = false;
		solver.tTimestep = 1000.0;
		solver.timeDelta = 1.0;
		
		solver.solve();
		
		
	}
}
