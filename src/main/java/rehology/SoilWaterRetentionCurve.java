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

import java.util.Map;

import physicalquantities.SoilParameters;
import physicalquantities.Variables;

/**
 * The SWRC abstract class.
 * @author Niccolo' Tubini
 *
 */

public abstract class SoilWaterRetentionCurve {
	
	private int ID;

	/**
	 * FIXME necessary??
	 * This method return the value of suction at which the derivative
	 * of moisture capacity with respect of suction has its maximum.
	 * @return
	 */
	public double getPsiStar1(int i){
		
		ID = SoilParameters.elementsLabel.get(i);

		return SoilParameters.psiStar1[ID];
	}
	
	
	
	/**
	 * This method compute the saturation degree given the suction value 
	 * @param i control volume ID
	 * @return
	 */
	public double saturationDegree(int i) {
		
		ID = SoilParameters.elementsLabel.get(i);

		return (waterContent(i)-SoilParameters.thetaR[ID])/(SoilParameters.thetaS[ID]-SoilParameters.thetaR[ID]);
	};
	
	
	
	
	/**
	 * This method compute the saturation degree given the suction value 
	 * @param i control volume ID
	 * @return
	 */
	public double saturationDegree(double suction, int i) {
		
		ID = SoilParameters.elementsLabel.get(i);

		return (waterContent(suction,i)-SoilParameters.thetaR[ID])/(SoilParameters.thetaS[ID]-SoilParameters.thetaR[ID]);
	};
	
	
	
	/**
	 * This method compute the water content given the suction value 
	 * @param i control volume ID
	 * @return
	 */
	public abstract double waterContent(int i);
	
	
	
	/**
	 * This method compute the water content given the suction value 
	 * @param suction
	 * @param i control volume ID
	 * @return
	 */
	public abstract double waterContent(double suction, int i);
	
	
		
	/**
	 * This method compute the value of the derivative of theta given the suction value
	 * @param i control volume ID
	 * @return
	 */
	public abstract double dWaterContent(int i);
	
	
	
	/**
	 * This method compute the value of the derivative of theta given the suction value
	 * @param suction
	 * @param i control volume ID
	 * @return
	 */
	public abstract double dWaterContent(double suction, int i);
	
	

	/**
	 * This method compute the integral of p function, Jordan decomposition (Casulli and Zanolli, 2010)
	 * @param suction
	 * @return
	 */
	public abstract double pIntegral(int i);
	
	
	/**
	 * This method compute the integral of p function, Jordan decomposition (Casulli and Zanolli, 2010)
	 * @param suction
	 * @return
	 */
	public abstract double pIntegral(double suction, int i);
	
	
	/**
	 * This method compute the integral of q function, Jordan decomposition (Casulli and Zanolli, 2010)
	 * @param suction
	 * @return
	 */
	public abstract double qIntegral(int i);
	
	
	
	/**
	 * This method compute the integral of q function, Jordan decomposition (Casulli and Zanolli, 2010)
	 * @param suction
	 * @return
	 */
	public abstract double qIntegral(double suction, int i);
	
	
	
	/**
	 * This method compute the p function, Jordan decomposition (Casulli and Zanolli, 2010)
	 * @param suction
	 * @return
	 */
	public abstract double p(int i);
	
	
	
	/**
	 * This method compute the p function, Jordan decomposition (Casulli and Zanolli, 2010)
	 * @param suction
	 * @return
	 */
	public abstract double p(double suction, int i);
	
	
	
	/**
	 * This method compute the q function, Jordan decomposition (Casulli and Zanolli, 2010)
	 * @param suction
	 * @return
	 */
	public abstract double q(int i);
	
	
	
	/**
	 * This method compute the p function, Jordan decomposition (Casulli and Zanolli, 2010)
	 * @param suction
	 * @return
	 */
	public abstract double q(double suction, int i);
	
	
}