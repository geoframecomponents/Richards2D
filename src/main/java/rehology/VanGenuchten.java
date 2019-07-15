/*
 * GNU GPL v3 License
 *
 * Copyright 2017  Niccolo` Tubini
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

import physicalquantities.SoilParameters;
import physicalquantities.Variables;

/**
 * Van Genuchten's SWRC model
 * @author Niccolo' Tubini
 */

public class VanGenuchten extends SoilWaterRetentionCurve {
	
	private double m;
	private int ID;

	
	
	//@Override
	public double waterContent(int i){
		
		ID = SoilParameters.elementsLabel.get(i);
		
		this.m = 1-1/SoilParameters.par1[ID];
		if(Variables.waterSuctions.get(i) <= 0) {
		    return SoilParameters.thetaR[ID] + (SoilParameters.thetaS[ID] - SoilParameters.thetaR[ID]) / Math.pow(1.0 + Math.pow(Math.abs(SoilParameters.par2[ID]*Variables.waterSuctions.get(i)), SoilParameters.par1[ID]), this.m);
		} else {
		    return  SoilParameters.thetaS[ID] + 9.81*( SoilParameters.alphaSpecificStorage[ID] + SoilParameters.thetaS[ID]*SoilParameters.betaSpecificStorage[ID])*Variables.waterSuctions.get(i) ;
		}

	}
	
	
	
	//@Override
	public double waterContent(double suction, int i){
		
		ID = SoilParameters.elementsLabel.get(i);
		
		this.m = 1-1/SoilParameters.par1[ID];
		if(suction <= 0) {
		    return SoilParameters.thetaR[ID] + (SoilParameters.thetaS[ID] - SoilParameters.thetaR[ID]) / Math.pow(1.0 + Math.pow(Math.abs(SoilParameters.par2[ID]*suction), SoilParameters.par1[ID]), this.m);
		} else {
		    return  SoilParameters.thetaS[ID] + 9.81*( SoilParameters.alphaSpecificStorage[ID] + SoilParameters.thetaS[ID]*SoilParameters.betaSpecificStorage[ID])*suction;
		}

	}
	
	
	
	//@Override
	public double dWaterContent(int i){
		
		ID = SoilParameters.elementsLabel.get(i);
		
		this.m = 1-1/SoilParameters.par1[ID];
		if (Variables.waterSuctions.get(i) <= 0) {
		    return SoilParameters.par2[ID]*SoilParameters.par1[ID]*this.m*(SoilParameters.thetaS[ID] - SoilParameters.thetaR[ID]) / Math.pow(1.0 + Math.pow(Math.abs(SoilParameters.par2[ID]*Variables.waterSuctions.get(i)), SoilParameters.par1[ID]), this.m + 1.0)*Math.pow(Math.abs(SoilParameters.par2[ID]*Variables.waterSuctions.get(i)), SoilParameters.par1[ID] - 1.0);
		} else {
		    return 9.81*( SoilParameters.alphaSpecificStorage[ID] + SoilParameters.thetaS[ID]*SoilParameters.betaSpecificStorage[ID] );
		}
		
	}
	
	
	//@Override
	public double dWaterContent(double suction, int i){
		
		ID = SoilParameters.elementsLabel.get(i);

		this.m = 1-1/SoilParameters.par1[ID];
		if (suction <= 0) {
		    return SoilParameters.par2[ID]*SoilParameters.par1[ID]*this.m*(SoilParameters.thetaS[ID] - SoilParameters.thetaR[ID]) / Math.pow(1.0 + Math.pow(Math.abs(SoilParameters.par2[ID]*suction), SoilParameters.par1[ID]), this.m + 1.0)*Math.pow(Math.abs(SoilParameters.par2[ID]*suction), SoilParameters.par1[ID] - 1.0);
		} else {
		    return 9.81*( SoilParameters.alphaSpecificStorage[ID] + SoilParameters.thetaS[ID]*SoilParameters.betaSpecificStorage[ID] );
		}
		
	}
	
	
	
	//@Override
	public double pIntegral(int i){
		
		ID = SoilParameters.elementsLabel.get(i);

		if(Variables.waterSuctions.get(i) <= SoilParameters.psiStar1[ID]) {
			return this.waterContent(i);
		} else if(SoilParameters.psiStar1[ID]<Variables.waterSuctions.get(i) && Variables.waterSuctions.get(i)<=0) {
			return this.waterContent(SoilParameters.psiStar1[ID],i) + this.dWaterContent(SoilParameters.psiStar1[ID],i)*(Variables.waterSuctions.get(i) - SoilParameters.psiStar1[ID]);
		} else {
			return  this.waterContent(SoilParameters.psiStar1[ID],i) + this.dWaterContent(SoilParameters.psiStar1[ID],i)*(Variables.waterSuctions.get(i) - SoilParameters.psiStar1[ID]) + this.dWaterContent(i)*Variables.waterSuctions.get(i);
		}

	}
	
	
	
	//@Override
	public double pIntegral(double suction, int i){
		
		ID = SoilParameters.elementsLabel.get(i);

		if(suction <= SoilParameters.psiStar1[ID]) {
			return this.waterContent(suction, i);
		} else if(SoilParameters.psiStar1[ID]<suction && suction<=0) {
			return this.waterContent(SoilParameters.psiStar1[ID],i) + this.dWaterContent(SoilParameters.psiStar1[ID],i)*(suction - SoilParameters.psiStar1[ID]);
		} else {
			return  this.waterContent(SoilParameters.psiStar1[ID],i) + this.dWaterContent(SoilParameters.psiStar1[ID],i)*(suction - SoilParameters.psiStar1[ID]) + this.dWaterContent(suction,i)*suction;
		}

	}
	
	
	
	//@Override
	public double qIntegral(int i){
		
		return pIntegral(i) - this.waterContent(i);

	}
	
	
	
	//@Override
	public double qIntegral(double suction, int i){
		
		return pIntegral(suction, i) - this.waterContent(suction, i);

	}
	
	
	
	//@Override
	public double p(int i){
		
		ID = SoilParameters.elementsLabel.get(i);

		if (Variables.waterSuctions.get(i) <= SoilParameters.psiStar1[ID]) {
			return this.dWaterContent(i);
		} else if (SoilParameters.psiStar1[ID]<Variables.waterSuctions.get(i) && Variables.waterSuctions.get(i)<=0) {
			return this.dWaterContent(SoilParameters.psiStar1[ID],i);
		} else {
			return this.dWaterContent(SoilParameters.psiStar1[ID],i) + this.dWaterContent(i);
		}

	}
	
	
	
	//@Override
	public double p(double suction, int i){
		
		ID = SoilParameters.elementsLabel.get(i);

		if (suction <= SoilParameters.psiStar1[ID]) {
			return this.dWaterContent(suction,i);
		} else if (SoilParameters.psiStar1[ID]<suction && suction<=0) {
			return this.dWaterContent(SoilParameters.psiStar1[ID],i);
		} else {
			return this.dWaterContent(SoilParameters.psiStar1[ID],i) + this.dWaterContent(suction,i);
		}

	}
	
	
	
	
	//@Override
	public double q(int i){
		
		return p(i) - this.dWaterContent(i);

	}
	
	
	
	//@Override
	public double q(double suction, int i){
		
		return p(suction,i) - this.dWaterContent(suction,i);

	}


}
