/*
 * GNU GPL v3 License
 *
 * Copyright 2019  Niccolo` Tubini
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

//import junit.framework.Assert;
import physicalquantities.*;
import rehology.*;
import topology.Topology;
/**
 * This class carries out the Nested-Newton algorithm
 * (A NESTED NEWTON-TYPE ALGORITHM FOR FINITE VOLUME METHODS SOLVING RICHARDS' EQUATION IN MIXED FORM, Casulli V., Zanolli P., Journal Scientific Computing, 2010)
 *  @author Niccolo' Tubini
 */

public class NestedNewton {
	private double outerResidual;
	private double innerResidual;

	int nestedNewton;
	int MAXITER_NEWT;
	//int NUM_CONTROL_VOLUMES;

	double newtonTolerance;
	double tmp;
	double timeDelta;
	//double[] psis;
	//double[] mainDiagonal;
	//double[] upperDiagonal;
	//double[] lowerDiagonal;
	//double[] rhss;
	//double[] dx;

	//double[] thetaS;
	//double[] thetaR;
	//double[] par1SWRC;
	//double[] par2SWRC;
	//double[] par3SWRC;
	//double[] par4SWRC;
	//double[] par5SWRC;
	//double[] alphaSpecificStorage;
	//double[] betaSpecificStorage;

	Map<Integer, Double> rhss;
	Map<Integer, Double> Apsi;
	Map<Integer, Double> fs;
	Map<Integer, Double> fks;
	//double[] bb;
	//double[] cc;
	Map<Integer, Double> dis;
	Map<Integer, Double> dpsis;
	Map<Integer, Double> psis_outer;
	Map<Integer, Double> psism;

	SoilWaterRetentionCurve swrc;
	//TotalDepth totalDepth;
	//Thomas thomasAlg = new Thomas();

	Matop matop;
	ConjugateGradientMethod cg;



	/**
	 * @param nestedNewton control parameter to choose between simple Newton method (0), or the nested Newton one (1)
	 * @param newtonTolerance prefixed tolerance representing the maximum mass balance error allowed  
	 * @param MAXITER_NEWT prefixed maximum number of iteration
	 * @param NUM_CONTROL_VOLUMES number of control volumes
	 * @param soilPar is the class to compute the soil hydraulic properties
	 * @param totalDepth is the class to compute the total water depth
	 * @param par1SWRC vector containing the first parameter of the SWRC, it is a vector of length NUM_CONTROL_VOLUMES-1
	 * @param par2SWRC vector containing the second parameter of the SWRC, it is a vector of length NUM_CONTROL_VOLUMES-1
	 * @param thetaR vector containing the adimensional residual water contentfor each control volume, it is a vector of length NUM_CONTROL_VOLUMES-1
	 * @param thetaS vector containing the adimensional water content at saturation for each control volume, it is a vector of length NUM_CONTROL_VOLUMES-1
	 */
	public NestedNewton(int nestedNewton, double newtonTolerance, int MAXITER_NEWT, SoilWaterRetentionCurve swrc, Matop matop, double cgTolerance){

		this.nestedNewton = nestedNewton;
		this.newtonTolerance = newtonTolerance;
		this.MAXITER_NEWT = MAXITER_NEWT;
		this.swrc = swrc;
		this.matop = matop;
		this.cg = new ConjugateGradientMethod(matop, cgTolerance);


		fs			  = new HashMap<Integer, Double>();
		fks			  = new HashMap<Integer, Double>();
		dis			  = new HashMap<Integer, Double>();
		dpsis		  = new HashMap<Integer, Double>();
		psis_outer	  = new HashMap<Integer, Double>();
		psism          = new HashMap<Integer, Double>();
	}



	/**
	 * @param psis vector contains the suction values, it is a vector of length NUM_CONTROL_VOLUMES
	 * @param upperDiagonal upper diagonal of the coefficient matrix A of the linear system, it is a vector of length NUM_CONTROL_VOLUMES
	 * @param mainDiagonal main diagonal of the coefficient matrix A of the linear system, it is a vector of length NUM_CONTROL_VOLUMES
	 * @param lowerDiagonal lower diagonal of the coefficient matrix A of the linear system, it is a vector of length NUM_CONTROL_VOLUMES
	 * @param rhss right hand side term of the linear system, it is a vector of length NUM_CONTROL_VOLUMES
	 */
	public void set(Map<Integer, Double> rhss){

		this.rhss = rhss;

	}



	public void solver(){



		// Initial guess for the outer iteration
		for(Integer element : Topology.s_i.keySet()) {
			//			if(i==NUM_CONTROL_VOLUMES-1) {
			//				psis[i] = Math.max(psis[i],0.1);
			//				//System.out.println(i +"   "+psis[i]);
			tmp = Math.min(Variables.waterSuctions.get(element), SoilParameters.psiStar1[SoilParameters.elementsLabel.get(element)]);
			Variables.waterSuctions.put(element, tmp);
			//System.out.println(i +"   "+psis[i]);
		}

		//// OUTER CYCLE ////
		for(int i = 0; i < MAXITER_NEWT; i++) {
			// I have to assign 0 to outerResidual otherwise I will take into account of the previous error
			outerResidual = 0.0;
			for(Integer element : Topology.s_i.keySet()) {
				dis.put(element, 0.0);
			}
			//Apsi = matop2DRichards(dis, Variables.waterSuctions);
			Apsi = matop.solve(dis, Variables.waterSuctions);
//			System.out.println("Apsi :");
//			for(Integer element : Topology.s_i.keySet()) {
//				System.out.println("\t"+ element + "\t" + Apsi.get(element));
//			}
//			System.out.println("fs :");
			for(Integer element : Topology.s_i.keySet()) {
				tmp =  swrc.dWaterContent(Variables.waterSuctions.get(element),element)*Geometry.elementsArea.get(element);
				dis.put(element, tmp );
				tmp = swrc.waterContent(Variables.waterSuctions.get(element),element)*Geometry.elementsArea.get(element) - rhss.get(element) + Apsi.get(element);
				fs.put(element, tmp);
//				System.out.println(element +": " +swrc.waterContent(Variables.waterSuctions.get(element),element) +" "+Geometry.elementsArea.get(element) +" "+rhss.get(element)+" "+ Apsi.get(element)+" " +fs.get(element) + " " + dis.get(element));
//				System.out.println("\t"+ element + "\t" + fs.get(element));
		
				outerResidual += tmp*tmp;
			}
			outerResidual = Math.pow(outerResidual,0.5);  
			System.out.println("\t\t-Outer iteration " + i + " with residual " +  outerResidual);
			if(outerResidual < newtonTolerance) {
				break;
			}
			if(nestedNewton == 0){

			}else{

				// Initial guess for the inner iteration (optional)
				for(Integer element : Topology.s_i.keySet()) {
					//			if(i==NUM_CONTROL_VOLUMES-1) {
					//				psis[i] = Math.max(psis[i],0.1);
					//				//System.out.println(i +"   "+psis[i]);
					psis_outer.put(element, Variables.waterSuctions.get(element));
					//Variables.waterSuctions.put(element, Math.max(Variables.waterSuctions.get(element), SoilParameters.psiStar1[SoilParameters.elementsLabel.get(element)]));
				}

				//// INNER CYCLE ////
				for(int j = 0; j < MAXITER_NEWT; j++) {
					// I have to assign 0 to innerResidual otherwise I will take into account of the previous error
					innerResidual = 0.0; 
					for(Integer element : Topology.s_i.keySet()) {
						dis.put(element, 0.0);
					}
					//Apsi = matop2DRichards(dis, Variables.waterSuctions);
					Apsi = matop.solve(dis, Variables.waterSuctions);
//					System.out.println("Apsi :");
//					for(Integer element : Topology.s_i.keySet()) {
//						System.out.println("\t"+ element + "\t" + Apsi.get(element));
//					}
//					System.out.println("fks :");
					for(Integer element : Topology.s_i.keySet()) {

						tmp = (swrc.p(Variables.waterSuctions.get(element),element) - swrc.q(psis_outer.get(element),element))*Geometry.elementsArea.get(element);
						dis.put(element, tmp);
						tmp = swrc.pIntegral(Variables.waterSuctions.get(element),element)*Geometry.elementsArea.get(element)
								- ( swrc.qIntegral(psis_outer.get(element),element) + swrc.q(psis_outer.get(element),element)*(Variables.waterSuctions.get(element) - psis_outer.get(element)) )*Geometry.elementsArea.get(element)
								- rhss.get(element) + Apsi.get(element);
						fks.put(element, tmp);

//						System.out.println(element +" "+fks.get(element) + " " + dis.get(element));
//						System.out.println("\t"+ element + "\t" + fks.get(element));

						innerResidual += tmp*tmp;
					}

					innerResidual = Math.pow(innerResidual,0.5);

					System.out.println("\t\t\t-Inner iteration " + j + " with residual " +  innerResidual);    

					if(innerResidual < newtonTolerance) {
//						System.out.println("Psi:");
//						for(Integer element : Topology.s_i.keySet()) {
//							System.out.println("\t"+ element + "\t" + Variables.waterSuctions.get(element));
//
//						}
//
//						System.out.println("\n\n");
//						System.out.println("dpsi:");
//						for(Integer element : Topology.s_i.keySet()) {
//							System.out.println("\t"+ element + "\t" + dpsis.get(element));
//
//						}
//						
						break;
					}
					//// CONJUGATE GRADIENT METHOD////
					//dpsis = conjugateGradientMethod(dis, fks);
					dpsis = cg.solve(dis, fks);
//					System.out.println("dpsis :");
//					for(Integer element : Topology.s_i.keySet()) {
//						System.out.println("\t"+ element + "\t" + dpsis.get(element));
//					}
					for(Integer element : Topology.s_i.keySet()) {
//						psism.put(element, Variables.waterSuctions.get(element));
						tmp  = Variables.waterSuctions.get(element)-dpsis.get(element);
						Variables.waterSuctions.put(element, tmp );
					}
//					if(j>1) {
//						for(Integer element : Topology.s_i.keySet()) {
//							Variables.waterSuctions.put(element, Math.min( Variables.waterSuctions.get(element), psism.get(element)) );
//							if(i>1) {
//									Variables.waterSuctions.put(element, Math.max( Variables.waterSuctions.get(element), psis_outer.get(element)) );
//							}
//						}
//					}
				} //// INNER CYCLE END ////

			}
			//		return psis;
//			if(i>1) {
//				for(Integer element : Topology.s_i.keySet()) {
//					Variables.waterSuctions.put(element, Math.max( Variables.waterSuctions.get(element), psis_outer.get(element)) );
//				}
//			}
		} //// OUTER CYCLE END ////

	}


	private Map<Integer, Double> matop2DRichards(Map<Integer, Double> dis, Map<Integer, Double> variable){

		Map<Integer, Double> Apsi = new HashMap<Integer,Double>();
		double sideFlux = 0.0;
		//double sigma;
		double kappa;

		/*
		 * Da cancellare meglio fare un ciclo sui lati e quindi calcolare i flussi.
		 * piu` veloce perche` lo stesso lato non lo devo calcolare due volte,
		 * non ho bisogno di calcolare la quantita sigma perche` quando vado ad aggiungere il flusso 
		 * all elemento decido io il segno in base a left e right ( vedi matop.f90 linee 92 e seguenti)
		 */
		//		for(Integer element : Topology.s_i.keySet()) {
		//			sideFlux = 0.0;
		//			for(Integer edge : Topology.s_i.get(element)) {
		//				/*
		//				 * FIXME: everywhere NO FLUX BOUNDARY CONDITION
		//				 */
		//				if(Topology.r.get(edge)==0) {
		//					sideFlux += 0.0; // depends on the type of the boundary condition and its value
		//				} else {
		//					sigma = (Topology.r.get(edge) - 2*element + Topology.l.get(edge))/(Topology.r.get(edge) - Topology.l.get(edge));
		//					kappa = Math.max( Variables.kappas.get(Topology.l.get(edge)), Variables.kappas.get(Topology.r.get(edge)) )*Geometry.edgesLenght.get(edge);
		//					sideFlux += kappa * (Variables.waterSuctions.get(Topology.l.get(edge))-Variables.waterSuctions.get(Topology.r.get(edge)))/Geometry.delta_j.get(edge);
		//				}
		//			}
		//			Apsi.put(element, sideFlux);
		//			
		//		}

		for(Integer element : Topology.s_i.keySet()) {
			Apsi.put(element, dis.get(element)*variable.get(element));
		}


		for(Integer edge : Topology.l.keySet()) {
			/*
			 * FIXME: everywhere NO FLUX BOUNDARY CONDITION
			 */
			sideFlux = 0.0;

			if(Topology.r.get(edge)==0) {
				sideFlux = 0.0;  // depends on the type of the boundary condition and its value
				Apsi.put(Topology.l.get(edge), Apsi.get(Topology.l.get(edge))-sideFlux);
			} else {
				kappa = Math.max( Variables.kappas.get(Topology.l.get(edge)), Variables.kappas.get(Topology.r.get(edge)) )*Geometry.edgesLenght.get(edge);
				sideFlux =  Variables.timeDelta * kappa * (variable.get(Topology.r.get(edge))-variable.get(Topology.l.get(edge)))/Geometry.delta_j.get(edge);
				Apsi.put( Topology.l.get(edge), Apsi.get(Topology.l.get(edge))-sideFlux );
				Apsi.put( Topology.r.get(edge), Apsi.get(Topology.r.get(edge))+sideFlux );
			}
		}
//		System.out.println("\n\nApsi:");
//		for(Integer element : Topology.s_i.keySet()) {
//			System.out.println("\t" + element + "\t" + Apsi.get(element));
//		}
		return Apsi;
	}


	private Map<Integer, Double> conjugateGradientMethod(Map<Integer, Double> dis, Map<Integer, Double> b){

		Map<Integer, Double> residual = new HashMap<Integer, Double>();
		Map<Integer, Double> x = new HashMap<Integer, Double>(b);
		Map<Integer, Double> p = new HashMap<Integer, Double>();
		Map<Integer, Double> Apsi = new HashMap<Integer, Double>();
		double alpha;
		double alphak;
		double lambda;
		double tmp;

		Apsi = matop2DRichards(dis, x);
		alpha = 0.0;
		for(Integer element : Topology.s_i.keySet()) {
			residual.put(element, x.get(element)-Apsi.get(element));
			p.put(element, residual.get(element));
			alpha += residual.get(element)*residual.get(element);
		}

		for(int k=1; k<4*x.size();k++) {

			if(Math.sqrt(alpha)<Math.pow(10, -12)) {
				break;
			}

			tmp = 0.0;
			Apsi = matop2DRichards(dis, p);
			for(Integer element : Topology.s_i.keySet()) {
				tmp += p.get(element)*Apsi.get(element);
			}
			lambda = alpha/tmp;

			alphak = alpha;
			alpha = 0.0;
			for(Integer element : Topology.s_i.keySet()) {
				x.put(element, x.get(element) + lambda*p.get(element));
				residual.put(element, residual.get(element) - lambda*Apsi.get(element));
				alpha += residual.get(element)*residual.get(element);
			}

			for(Integer element : Topology.s_i.keySet()) {
				p.put(element, residual.get(element)+alpha/alphak*p.get(element) );
			}

		}

		return x;
	}

}
