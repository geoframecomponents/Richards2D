package rehology;

import swrc.SoilWaterRetentionCurve;

public abstract class UnsaturatedHydraulicConductivity {

	protected double saturationDegree;
	protected SoilWaterRetentionCurve modelSWRC;
	
	
	
	UnsaturatedHydraulicConductivity(SoilWaterRetentionCurve modelSWRC) {
		this.modelSWRC = modelSWRC;
	}
	
	
	
	/**
	 * This method compute the hydraulic conductivity.
	 * @param suction
	 * @return
	 */
	public abstract double hydraulicConductivity(int i);
	
	
	
	/**
	 * This method compute the hydraulic conductivity.
	 * @param suction
	 * @return
	 */
	public abstract double hydraulicConductivity(double suction, int i);
	
	
	
	/**
	 * This method compute the derivative of hydraulic conductivity with respect to 
	 * the water content.
	 * @param suction
	 * @return
	 */
	public abstract double dHydraulicConductivity(double suction, int i);

}
