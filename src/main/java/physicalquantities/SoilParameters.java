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

package physicalquantities;

import java.util.HashMap;
import java.util.Map;

public class SoilParameters {
	
	static SoilParameters uniqueInstance;

	public static SoilParameters getInstance() {
		/*if (uniqueInstance == null) {
				uniqueInstance = new Variables(waterSuction, temperature);
			}*/
		return uniqueInstance;
	}

	public static SoilParameters getInstance(Map<Integer, Integer> elementsLabel, double[] alphaSpecificStorage, double[] betaSpecificStorage, double[] kappaSaturation, 
			double[] par1, double[] par2, double[] par3, double[] par4, double[] par5, 
			double[] psiStar1, double[] psiStar2, double[] psiStar3, double[] thetaR, double[] thetaS) {
		if (uniqueInstance == null) {
			uniqueInstance = new SoilParameters ( elementsLabel, alphaSpecificStorage, betaSpecificStorage, kappaSaturation, 
					par1, par2, par3, par4,  par5, 
					psiStar1, psiStar2, psiStar3, thetaR, thetaS);
		}
		return uniqueInstance;
	}

	public static Map<Integer, Integer> elementsLabel;
	public static double[] alphaSpecificStorage;
	public static double[] betaSpecificStorage;
	public static double[] kappaSaturation;
	public static double[] par1;
	public static double[] par2;
	public static double[] par3;
	public static double[] par4;
	public static double[] par5;
	public static double[] psiStar1;
	public static double[] psiStar2;
	public static double[] psiStar3;
	public static double[] thetaR;
	public static double[] thetaS;


	public SoilParameters(Map<Integer, Integer> elementsLabel, double[] alphaSpecificStorage, double[] betaSpecificStorage, double[] kappaSaturation, 
			double[] par1, double[] par2, double[] par3, double[] par4, double[] par5, 
			double[] psiStar1, double[] psiStar2, double[] psiStar3, double[] thetaR, double[] thetaS) {

		SoilParameters.elementsLabel = new HashMap<Integer, Integer>();
		SoilParameters.elementsLabel = elementsLabel; 
		SoilParameters.alphaSpecificStorage = alphaSpecificStorage;
		SoilParameters.betaSpecificStorage = betaSpecificStorage;
		SoilParameters.kappaSaturation = kappaSaturation;
		SoilParameters.par1 = par1;
		SoilParameters.par2 = par2;
		SoilParameters.par3 = par3;
		SoilParameters.par4 = par4;
		SoilParameters.par5 = par5;
		SoilParameters.psiStar1 = psiStar1;
		SoilParameters.psiStar2 = psiStar2;
		SoilParameters.psiStar3 = psiStar3;
		SoilParameters.thetaR = thetaR;
		SoilParameters.thetaS = thetaS;

	}


}
