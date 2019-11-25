package edu.scripps.yates.nucleome.turboID;

public enum TurboID_Channel_Norm {
	nmT2_A1(Replicate.A, TurboIDExperimentType.EMD), //
	nmT3_A2(Replicate.A, TurboIDExperimentType.EMD), //
	nmT2_A7(Replicate.Arerun, TurboIDExperimentType.EMD), //
	nmT3_A8(Replicate.Arerun, TurboIDExperimentType.EMD), //
	nmT2_A3(Replicate.B, TurboIDExperimentType.EMD), //
	nmT3_A4(Replicate.B, TurboIDExperimentType.EMD), //
	nmT4_B1(Replicate.A, TurboIDExperimentType.LBR), //
	nmT5_B2(Replicate.A, TurboIDExperimentType.LBR), //
	nmT4_B7(Replicate.Arerun, TurboIDExperimentType.LBR), //
	nmT5_B8(Replicate.Arerun, TurboIDExperimentType.LBR), //
	nmT4_B3(Replicate.B, TurboIDExperimentType.LBR), //
	nmT5_B4(Replicate.B, TurboIDExperimentType.LBR), //
	nmT6_C1(Replicate.A, TurboIDExperimentType.SUN1), //
	nmT7_C2(Replicate.A, TurboIDExperimentType.SUN1), //
	nmT6_C7(Replicate.Arerun, TurboIDExperimentType.SUN1), //
	nmT7_C8(Replicate.Arerun, TurboIDExperimentType.SUN1), //
	nmT6_C3(Replicate.B, TurboIDExperimentType.SUN1), //
	nmT7_C4(Replicate.B, TurboIDExperimentType.SUN1), //
	nmT8_D1(Replicate.A, TurboIDExperimentType.MAN1), //
	nmT9_D2(Replicate.A, TurboIDExperimentType.MAN1), //
	nmT8_D7(Replicate.Arerun, TurboIDExperimentType.MAN1), //
	nmT9_D8(Replicate.Arerun, TurboIDExperimentType.MAN1), //
	nmT8_D3(Replicate.B, TurboIDExperimentType.MAN1), //
	nmT9_D4(Replicate.B, TurboIDExperimentType.MAN1), //
	nmT10_E1(Replicate.A, TurboIDExperimentType.TURBOID_ONLY), //
	nmT11_E2(Replicate.A, TurboIDExperimentType.TURBOID_ONLY), //
	nmT10_E7(Replicate.Arerun, TurboIDExperimentType.TURBOID_ONLY), //
	nmT11_E8(Replicate.Arerun, TurboIDExperimentType.TURBOID_ONLY), //
	nmT10_E3(Replicate.B, TurboIDExperimentType.TURBOID_ONLY), //
	nmT11_E4(Replicate.B, TurboIDExperimentType.TURBOID_ONLY), //
	nmT1_mix_a(Replicate.A, TurboIDExperimentType.MIX), //
	nmT1_mix_arerun(Replicate.Arerun, TurboIDExperimentType.MIX), //
	nmT1_mix_b(Replicate.B, TurboIDExperimentType.MIX), //
	nmT1_mix_brerun(Replicate.Brerun, TurboIDExperimentType.MIX), //
	nmT2_A5(Replicate.Brerun, TurboIDExperimentType.EMD), //
	nmT3_A6(Replicate.Brerun, TurboIDExperimentType.EMD), //
	nmT4_B5(Replicate.Brerun, TurboIDExperimentType.LBR), //
	nmT5_B6(Replicate.Brerun, TurboIDExperimentType.LBR), //
	nmT6_C5(Replicate.Brerun, TurboIDExperimentType.SUN1), //
	nmT7_C6(Replicate.Brerun, TurboIDExperimentType.SUN1), //
	nmT8_D5(Replicate.Brerun, TurboIDExperimentType.MAN1), //
	nmT9_D6(Replicate.Brerun, TurboIDExperimentType.MAN1), //
	nmT10_E5(Replicate.Brerun, TurboIDExperimentType.TURBOID_ONLY), //
	nmT11_E6(Replicate.Brerun, TurboIDExperimentType.TURBOID_ONLY), //

	norm_Cy_A1(Replicate.CyA, TurboIDExperimentType.EMD), //
	norm_Cy_A2(Replicate.CyA, TurboIDExperimentType.EMD), //
	norm_Cy_A3(Replicate.CyB, TurboIDExperimentType.EMD), //
	norm_Cy_A4(Replicate.CyB, TurboIDExperimentType.EMD), //
	norm_Cy_B1(Replicate.CyA, TurboIDExperimentType.LBR), //
	norm_Cy_B2(Replicate.CyA, TurboIDExperimentType.LBR), //
	norm_Cy_B3(Replicate.CyB, TurboIDExperimentType.LBR), //
	norm_Cy_B4(Replicate.CyB, TurboIDExperimentType.LBR), //
	norm_Cy_C1(Replicate.CyA, TurboIDExperimentType.SUN1), //
	norm_Cy_C2(Replicate.CyA, TurboIDExperimentType.SUN1), //
	norm_Cy_C3(Replicate.CyB, TurboIDExperimentType.SUN1), //
	norm_Cy_C4(Replicate.CyB, TurboIDExperimentType.SUN1), //
	norm_Cy_D1(Replicate.CyA, TurboIDExperimentType.MAN1), //
	norm_Cy_D2(Replicate.CyA, TurboIDExperimentType.MAN1), //
	norm_Cy_D3(Replicate.CyB, TurboIDExperimentType.MAN1), //
	norm_Cy_D4(Replicate.CyB, TurboIDExperimentType.MAN1), //
	norm_Cy_E1(Replicate.CyA, TurboIDExperimentType.TURBOID_ONLY), //
	norm_Cy_E2(Replicate.CyA, TurboIDExperimentType.TURBOID_ONLY), //
	norm_Cy_E3(Replicate.CyB, TurboIDExperimentType.TURBOID_ONLY), //
	norm_Cy_E4(Replicate.CyB, TurboIDExperimentType.TURBOID_ONLY), //
	norm_CyRa(Replicate.CyA, TurboIDExperimentType.MIX), //
	norm_CyRb(Replicate.CyB, TurboIDExperimentType.MIX) //
	;
	private final Replicate replicate;
	private final TurboIDExperimentType expType;

	private TurboID_Channel_Norm(Replicate replicate, TurboIDExperimentType expType) {
		this.replicate = replicate;
		this.expType = expType;
	}

	public Replicate getReplicate() {
		return replicate;
	}

	public static TurboID_Channel_Norm getByName(String name) {
		for (final TurboID_Channel_Norm value : values()) {
			if (value.name().equals(name)) {
				return value;
			}
		}
		return null;
	}

	public TurboIDExperimentType getExpType() {
		return expType;
	}
}
