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

import java.net.URISyntaxException;
import java.util.HashMap;
import java.util.Map;

import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorReader;
import org.junit.Test;

import bidimensionalProblemTimeDependent.WriteNetCDFRichards2D;
import bufferWriter.RichardsBuffer2D;
import generatemesh.GenerateMesh;
//import meshtopology.TopologyTriangularMesh;
//import monodimensionalProblemTimeDependent.WriteNetCDFRichards1D;
import readtriangularization.*;

public class TestCallRichardsSolver {

	@Test
	public void Test() throws Exception {

		//String fileName = "resources/input/square_with_subdomain_100.msh";
		//String fileName = "resources/input/square22_1_simmetrico.msh";
		//String fileName = "resources/input/layer_DD.msh";
		//String fileName = "resources/input/CasulliOK.csv";//100x60_NN_1.csv"; //column160.csv";
		String fileName = "resources/input/square23_2.msh";
		String splitter = " ";

		Readmsh reader = new Readmsh();
		//Readcsv reader = new Readcsv();
		reader.fileName = fileName;
		reader.splitter = splitter;
		reader.checkData = false;
		reader.process();
		
		//GenerateMeshTriangles generateMesh = new GenerateMeshTriangles();
		GenerateMesh generateMesh = new GenerateMesh();
		generateMesh.verticesCoordinates = reader.verticesCoordinates;
		generateMesh.elementsVertices = reader.elementsVertices;
		generateMesh.borderEdgesVertices = reader.borderEdgesVertices;
		generateMesh.borderEdgesLabel = reader.borderEdgesLabel;
		generateMesh.checkData = false;
		generateMesh.meshType = "triangular";
		generateMesh.geometryType = "EuclideanCartesian";
		generateMesh.process();
		
		//Boundary conditions
		String startDate = "2017-01-01 00:00";
		String endDate = "2017-01-01 00:05";
		int timeStepMinutes = 5;
		String fId = "ID";
		
		//NetCDF
		//String pathOutput = "C:\\Users\\Niccolo\\eclipse-workspace\\Richards2D\\resources\\MCBride2016.nc";
		String pathOutput = "C:\\Users\\Niccolo\\eclipse-workspace\\Richards2D\\resources\\square23_2.nc";
		String outputDescription = "";//Mc Bride 2016 problem 4. picard iteration 2, Nested Newton tol 1E-9, CG tol 1E-9, delta t 3600s.";
		
		
		String path ="resources/input/Test.csv";
				
		
		OmsTimeSeriesIteratorReader readerBC = getTimeseriesReader(path, fId, startDate, endDate, timeStepMinutes);
		
		RichardsBuffer2D buffer = new RichardsBuffer2D();
		WriteNetCDFRichards2D writeNetCDF = new WriteNetCDFRichards2D();
		
		double[] alphaSpecificStorage = new double[] {0.0, 0.0};
		double[] betaSpecificStorage = new double[] {0.0, 0.0};
		double[] ks = new double[] {0.0015167, 0.0000626};//new double[] {0.0000015167, 0.0000626};//0.00626};//{0.0000028, 0.0000123};0.00015167
		double[] par1SWRC = new double[] {1.3954, 2.239}; //1.76
		double[] par2SWRC = new double[] {1.04, 2.8};   //2.6
		double[] par3SWRC = null;
		double[] par4SWRC = null;
		double[] par5SWRC = null;
		double[] psiStar1 = new double[] {-1/1.04*Math.pow(( 0.3954/1.3954),1/1.3954), -1/2.8*Math.pow(( 1.239/2.239),1/2.239)}; //-0.431-0.274199255289329
		double[] psiStar2 = null;
		double[] psiStar3 = null;
		double[] thetaS = new double[] {0.4686, 0.3658};//0.4686};
		double[] thetaR = new double[] {0.1060, 0.0286};//0.1060};
		//double[] thetaR = new double[] {0.2262, 0.07818};//0.1060}; casulli2
		//double[] thetaR = new double[] {0.1060, 0.07818};//0.1060}; casulli4


		
		CallRichards2DSolver solver = new CallRichards2DSolver();
		solver.l = generateMesh.l;
		solver.r = generateMesh.r;
		solver.s_i = generateMesh.s_i;		
		solver.elementsArea = generateMesh.elementsArea;
		solver.edgesLenght = generateMesh.edgesLength;
		solver.delta_j = generateMesh.delta_j;
		solver.edgeNormalVector = generateMesh.edgeNormalVector;
		solver.elementsCentroidsCoordinates = generateMesh.elementsCentroidsCoordinates;
		solver.edgesCentroidsCoordinates = generateMesh.edgesCentroidsCoordinates;
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
		solver.edgesBoundaryBCType = generateMesh.edgeBoundaryBCType;
		solver.edgesBoundaryBCValue = generateMesh.edgeBoundaryBCValue;
		solver.soilHydraulicModel = "VanGenuchten";
		solver.typeUHCModel = "MualemVanGenuchten";
		solver.typeMatop = "2DRichards";
		solver.initialConditionType = "hydrostatic";
		solver.newtonTolerance = 0.000000001;
		solver.MAXITER_NEWT = 30;
		solver.picardIteration = 1; 
		solver.cgTolerance = 0.000000001;
		solver.checkData = false;
		solver.tTimestep = 300.0;
		solver.timeDelta = 300.0;
		solver.zSeePage = -2.0;
		
		while( readerBC.doProcess  ) {
			
			readerBC.nextRecord();	
			HashMap<Integer, double[]> bCValueMap = readerBC.outData;
			solver.inBC= bCValueMap;
			
			solver.inCurrentDate = readerBC.tCurrent;
			
			solver.solve();
			
			buffer.inputDate = readerBC.tCurrent;
			buffer.inputSpatialCoordinate = generateMesh.elementsCentroidsCoordinates;
			buffer.inputDualSpatialCoordinate = generateMesh.edgesCentroidsCoordinates;
			buffer.inputVariable = solver.outputToBuffer;
			
			buffer.solve();

			writeNetCDF.fileName = pathOutput;
			writeNetCDF.briefDescritpion = outputDescription;
			writeNetCDF.myVariables = buffer.myVariable;
			writeNetCDF.mySpatialCoordinateX = buffer.mySpatialCoordinateX;
			writeNetCDF.mySpatialCoordinateZ = buffer.mySpatialCoordinateZ;
			writeNetCDF.myDualSpatialCoordinateX = buffer.myDualSpatialCoordinateX;	
			writeNetCDF.myDualSpatialCoordinateZ = buffer.myDualSpatialCoordinateZ;	
			writeNetCDF.doProcess = readerBC.doProcess;
			writeNetCDF.writeNetCDF();
			
		}
		
		
	}
	
	private OmsTimeSeriesIteratorReader getTimeseriesReader( String inPath, String id, String startDate, String endDate,
			int timeStepMinutes ) throws URISyntaxException {
		OmsTimeSeriesIteratorReader reader = new OmsTimeSeriesIteratorReader();
		reader.file = inPath;
		reader.idfield = "ID";
		reader.tStart = startDate;
		reader.tTimestep = timeStepMinutes;
		reader.tEnd = endDate;
		reader.fileNovalue = "-9999";
		reader.initProcess();
		return reader;
	}
	
}
