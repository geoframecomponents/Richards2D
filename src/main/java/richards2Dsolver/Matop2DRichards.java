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

import physicalquantities.Geometry;
import physicalquantities.Variables;
import topology.Topology;

public class Matop2DRichards extends Matop {
	
	private static Matop2DRichards uniqueInstance;
	public Map<Integer, Double> Apsi;
	
	public static Matop2DRichards getInstance() {
		if (uniqueInstance == null) {
			uniqueInstance = new Matop2DRichards();
		}
		return uniqueInstance;
	}
	
	public Matop2DRichards() {
		Apsi = new HashMap<Integer,Double>();
	}

	public Map<Integer, Double> solve(Map<Integer, Double> dis, Map<Integer, Double> variable) {

		//Map<Integer, Double> Apsi = new HashMap<Integer,Double>();
		double sideFlux = 0.0;
		double sigma;
		double kappa;

		Apsi.clear();
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
				//sideFlux = Variables.timeDelta * kappa * (variable.get(Topology.r.get(edge))-variable.get(Topology.l.get(edge)))/Geometry.delta_j.get(edge);
				sideFlux = Variables.timeDelta * kappa * ( variable.get(Topology.r.get(edge))-variable.get(Topology.l.get(edge)) )/Geometry.delta_j.get(edge);
				Apsi.put( Topology.l.get(edge), Apsi.get(Topology.l.get(edge))-sideFlux );
				Apsi.put( Topology.r.get(edge), Apsi.get(Topology.r.get(edge))+sideFlux );
			}
		}
		//			System.out.println("\n\nApsi:");
		//			for(Integer element : Topology.s_i.keySet()) {
		//				System.out.println("\t" + element + "\t" + Apsi.get(element));
		//			}
		
		for(Integer edge : Topology.edgesBoundaryBCType.keySet()) {
			/*
			 * FIXME: boundary conditions
			 * if 0 no flux
			 * if 1 neumann 
			 * if 2 dirichlet
			 */
			
			if(Topology.edgesBoundaryBCType.get(edge)==2) { // dirichlet
				kappa = Variables.kappas.get(Topology.l.get(edge))*Geometry.edgesLenght.get(edge);
				sideFlux = Variables.timeDelta * kappa * ( -variable.get(Topology.l.get(edge)) )/Geometry.delta_j.get(edge);
				
				Apsi.put( Topology.l.get(edge), Apsi.get(Topology.l.get(edge))-sideFlux );
			}
			
		}

		return Apsi;
		
	}
	
}
