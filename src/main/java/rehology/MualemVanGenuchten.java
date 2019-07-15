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

package rehology;

import physicalquantities.Variables;
import physicalquantities.SoilParameters;

public class MualemVanGenuchten extends UnsaturatedHydraulicConductivity {
	
	private double m = -999.0;
	private int ID;
	
	
	MualemVanGenuchten(SoilWaterRetentionCurve modelSWRC) {
		super(modelSWRC);
		// TODO Auto-generated constructor stub
	}

	
	
	@Override
	public double hydraulicConductivity(int i) {
		
		ID = SoilParameters.elementsLabel.get(i);
		
		this.m = 1-1/SoilParameters.par1[ID];
		//super.saturationDegree = super.modelSWRC.saturationDegree(i); 
		super.saturationDegree = (Variables.thetasOld.get(i)-SoilParameters.thetaR[ID])/(SoilParameters.thetaS[ID]-SoilParameters.thetaR[ID]); 
		if(super.saturationDegree<1) {
			return SoilParameters.kappaSaturation[ID] * Math.pow(super.saturationDegree, 0.5 ) * Math.pow(1.0 - Math.pow(1.0 - Math.pow(super.saturationDegree, 1.0/this.m), this.m), 2.0);
		} else {
			return SoilParameters.kappaSaturation[ID];
		}

	}

	
	
	@Override
	public double hydraulicConductivity(double suction, int i) {
		
		ID = SoilParameters.elementsLabel.get(i);

		this.m = 1-1/SoilParameters.par1[ID];
		super.saturationDegree = super.modelSWRC.saturationDegree(suction, i); 
		if(super.saturationDegree<1) {
			return SoilParameters.kappaSaturation[ID] * Math.pow(super.saturationDegree, 0.5 ) * Math.pow(1.0 - Math.pow(1.0 - Math.pow(super.saturationDegree, 1.0/this.m), this.m), 2.0);
		} else {
			return SoilParameters.kappaSaturation[ID];
		}

	}

	
	
	@Override
	public double dHydraulicConductivity(double suction, int i) {
		
		ID = SoilParameters.elementsLabel.get(i);

		this.m = 1-1/SoilParameters.par1[ID];
		super.saturationDegree = super.modelSWRC.saturationDegree(suction, i); 
		if(super.saturationDegree<1) {
			return 0.5*SoilParameters.kappaSaturation[ID]/(SoilParameters.thetaS[ID]-SoilParameters.thetaR[ID]) * Math.pow( (super.saturationDegree), -0.5 ) * Math.pow( 1.0 - Math.pow( 1.0 - Math.pow( super.saturationDegree, 1.0/this.m ), this.m ), 2.0 ) + 
					SoilParameters.kappaSaturation[ID]*Math.pow( (super.saturationDegree), 0.5 ) * 2 * ( 1.0 - Math.pow( 1.0 - Math.pow( super.saturationDegree, 1.0/this.m ), this.m ) ) * ( this.m* Math.pow( (1-Math.pow( super.saturationDegree,1.0/this.m) ) ,this.m-1.0) ) * (1/this.m*Math.pow(super.saturationDegree, 1.0/this.m-1)); 
		} else {
			return 0.5*SoilParameters.kappaSaturation[ID]/(SoilParameters.thetaS[ID]-SoilParameters.thetaR[ID]);
		}

	}

}
