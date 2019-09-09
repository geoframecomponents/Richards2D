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

import topology.Topology;

/**
 * <h1>Conjugate gradient method </h1>
 * Matrix-free conjugate method for the solution of A*x=b
 * A is symmetric and positive definite.
 * <p>
 * 
 * We do not need to store the matrix, just to provide a class which computes
 * the matrix-vector product A*x
 * This implementation works for the 1D, 2D as well as 3D case. What changes, is the 
 * class that computes the matrix-vector.
 * 
 * @author Niccolo' Tubini, Michael Dumbser, and Riccardo Rigon
 * @version 0.1
 * @since 2018-11-23
 * @see <a href="https://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf">J. R. Shewchuk, An Introduction to the
 *	Conjugate Gradient Method Without the Agonizing Pain</a>
 */

public class ConjugateGradientMethod {
	
	double alpha;
	double alphak;
	double lambda;
	double tmp;
	double cgTolerance;
	Map<Integer, Double> residual;
	Map<Integer, Double> x;
	Map<Integer, Double> p;
	Map<Integer, Double> Apsi;
	Matop matop;
	
	public ConjugateGradientMethod(Matop matop, double cgTolerance) {
		
		residual = new HashMap<Integer, Double>();
		x = new HashMap<Integer, Double>();
		p = new HashMap<Integer, Double>();
		Apsi = new HashMap<Integer, Double>();
		this.matop = matop;
		this.cgTolerance = cgTolerance;
	}

	public Map<Integer, Double> solve(Map<Integer, Double> dis, Map<Integer, Double> b){

		x = b;

		Apsi = matop.solve(dis, x);
		alpha = 0.0;
		for(Integer element : Topology.s_i.keySet()) {
			residual.put(element, b.get(element)-Apsi.get(element));
			p.put(element, residual.get(element));
			alpha += residual.get(element)*residual.get(element);
		}
		int iter =1;
		for(int k=1; k<4*x.size();k++) {

			if(Math.sqrt(alpha)<=cgTolerance) {
//				System.out.println("\t\t\t\t\tk: " + k +"\tsqrt(alpha): " +Math.sqrt(alpha));
				break;
			}
//			System.out.println("\t\tk: " + k);
//			if(k==17) {
//				System.out.println("\t\tk = 17");
//			}
			tmp = 0.0;

			Apsi = matop.solve(dis, p);
//			if(k>=1700000) {
//				System.out.println("Apsi :");
//				for(Integer element : Topology.s_i.keySet()) {
//					System.out.println("\t"+ element + "\t" + Apsi.get(element));
//				}
//			}
			tmp = 0.0;
			for(Integer element : Topology.s_i.keySet()) {
				tmp += p.get(element)*Apsi.get(element);
			}
			lambda = alpha/tmp;
//			System.out.println("\t\t\t\ttmp: " + tmp + "\tlambda: " + lambda);

			alphak = alpha;
			alpha = 0.0;
			for(Integer element : Topology.s_i.keySet()) {
				tmp = x.get(element) + lambda*p.get(element);
				x.put(element, tmp);
				tmp = residual.get(element) - lambda*Apsi.get(element);
				residual.put(element, tmp);
				//alpha += residual.get(element)*residual.get(element);
				alpha += tmp*tmp;
			}
			
//			if(k>=170000) {
//				System.out.println("\n\nx :");
//				for(Integer element : Topology.s_i.keySet()) {
//					System.out.println("\t"+ element + "\t" + x.get(element));
//				}
//			}
//			if(k>=170000) {
//				System.out.println("\n\nresidual :");
//				for(Integer element : Topology.s_i.keySet()) {
//					System.out.println("\t"+ element + "\t" + residual.get(element));
//				}
//			}
//			System.out.println("\t\t\t\t\talpha: " + alpha + "\talphak: " + alphak);
			for(Integer element : Topology.s_i.keySet()) {
				tmp = residual.get(element)+alpha/alphak*p.get(element);
				p.put(element, tmp );
			}
//			if(k>=170000) {
//				System.out.println("\n\np :");
//				for(Integer element : Topology.s_i.keySet()) {
//					System.out.println("\t"+ element + "\t" + p.get(element));
//				}
//			}
			iter++;
		}
		System.out.println("\t\t\t\t\tk: " + iter +"\tsqrt(alpha): " +Math.sqrt(alpha));

		return x;
	}

}
