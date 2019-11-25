package edu.scripps.yates.nucleome.turboID;

public enum TurboID_Channel_Ori {
	T2_A1(Replicate.A, TurboIDExperimentType.EMD), //
	T3_A2(Replicate.A, TurboIDExperimentType.EMD), //
	T2_A7(Replicate.Arerun, TurboIDExperimentType.EMD), //
	T3_A8(Replicate.Arerun, TurboIDExperimentType.EMD), //
	T2_A3(Replicate.B, TurboIDExperimentType.EMD), //
	T3_A4(Replicate.B, TurboIDExperimentType.EMD), //
	T4_B1(Replicate.A, TurboIDExperimentType.LBR), //
	T5_B2(Replicate.A, TurboIDExperimentType.LBR), //
	T4_B7(Replicate.Arerun, TurboIDExperimentType.LBR), //
	T5_B8(Replicate.Arerun, TurboIDExperimentType.LBR), //
	T4_B3(Replicate.B, TurboIDExperimentType.LBR), //
	T5_B4(Replicate.B, TurboIDExperimentType.LBR), //
	T6_C1(Replicate.A, TurboIDExperimentType.SUN1), //
	T7_C2(Replicate.A, TurboIDExperimentType.SUN1), //
	T6_C7(Replicate.Arerun, TurboIDExperimentType.SUN1), //
	T7_C8(Replicate.Arerun, TurboIDExperimentType.SUN1), //
	T6_C3(Replicate.B, TurboIDExperimentType.SUN1), //
	T7_C4(Replicate.B, TurboIDExperimentType.SUN1), //
	T8_D1(Replicate.A, TurboIDExperimentType.MAN1), //
	T9_D2(Replicate.A, TurboIDExperimentType.MAN1), //
	T8_D7(Replicate.Arerun, TurboIDExperimentType.MAN1), //
	T9_D8(Replicate.Arerun, TurboIDExperimentType.MAN1), //
	T8_D3(Replicate.B, TurboIDExperimentType.MAN1), //
	T9_D4(Replicate.B, TurboIDExperimentType.MAN1), //
	T10_E1(Replicate.A, TurboIDExperimentType.TURBOID_ONLY), //
	T11_E2(Replicate.A, TurboIDExperimentType.TURBOID_ONLY), //
	T10_E7(Replicate.Arerun, TurboIDExperimentType.TURBOID_ONLY), //
	T11_E8(Replicate.Arerun, TurboIDExperimentType.TURBOID_ONLY), //
	T10_E3(Replicate.B, TurboIDExperimentType.TURBOID_ONLY), //
	T11_E4(Replicate.B, TurboIDExperimentType.TURBOID_ONLY), //
	T1_mix_a(Replicate.A, TurboIDExperimentType.MIX), //
	T1_mix_arerun(Replicate.Arerun, TurboIDExperimentType.MIX), //
	T1_mix_b(Replicate.B, TurboIDExperimentType.MIX), //
	T1_mix_brerun(Replicate.Brerun, TurboIDExperimentType.MIX), //
	T2_A5(Replicate.Brerun, TurboIDExperimentType.EMD), //
	T3_A6(Replicate.Brerun, TurboIDExperimentType.EMD), //
	T4_B5(Replicate.Brerun, TurboIDExperimentType.LBR), //
	T5_B6(Replicate.Brerun, TurboIDExperimentType.LBR), //
	T6_C5(Replicate.Brerun, TurboIDExperimentType.SUN1), //
	T7_C6(Replicate.Brerun, TurboIDExperimentType.SUN1), //
	T8_D5(Replicate.Brerun, TurboIDExperimentType.MAN1), //
	T9_D6(Replicate.Brerun, TurboIDExperimentType.MAN1), //
	T10_E5(Replicate.Brerun, TurboIDExperimentType.TURBOID_ONLY), //
	T11_E6(Replicate.Brerun, TurboIDExperimentType.TURBOID_ONLY), //

	ori_Cy_A1(Replicate.CyA, TurboIDExperimentType.EMD), //
	ori_Cy_A2(Replicate.CyA, TurboIDExperimentType.EMD), //
	ori_Cy_A3(Replicate.CyB, TurboIDExperimentType.EMD), //
	ori_Cy_A4(Replicate.CyB, TurboIDExperimentType.EMD), //
	ori_Cy_B1(Replicate.CyA, TurboIDExperimentType.LBR), //
	ori_Cy_B2(Replicate.CyA, TurboIDExperimentType.LBR), //
	ori_Cy_B3(Replicate.CyB, TurboIDExperimentType.LBR), //
	ori_Cy_B4(Replicate.CyB, TurboIDExperimentType.LBR), //
	ori_Cy_C1(Replicate.CyA, TurboIDExperimentType.SUN1), //
	ori_Cy_C2(Replicate.CyA, TurboIDExperimentType.SUN1), //
	ori_Cy_C3(Replicate.CyB, TurboIDExperimentType.SUN1), //
	ori_Cy_C4(Replicate.CyB, TurboIDExperimentType.SUN1), //
	ori_Cy_D1(Replicate.CyA, TurboIDExperimentType.MAN1), //
	ori_Cy_D2(Replicate.CyA, TurboIDExperimentType.MAN1), //
	ori_Cy_D3(Replicate.CyB, TurboIDExperimentType.MAN1), //
	ori_Cy_D4(Replicate.CyB, TurboIDExperimentType.MAN1), //
	ori_Cy_E1(Replicate.CyA, TurboIDExperimentType.TURBOID_ONLY), //
	ori_Cy_E2(Replicate.CyA, TurboIDExperimentType.TURBOID_ONLY), //
	ori_Cy_E3(Replicate.CyB, TurboIDExperimentType.TURBOID_ONLY), //
	ori_Cy_E4(Replicate.CyB, TurboIDExperimentType.TURBOID_ONLY), //
	ori_CyRa(Replicate.CyA, TurboIDExperimentType.MIX), //
	ori_CyRb(Replicate.CyB, TurboIDExperimentType.MIX) //

	;
	private final Replicate replicate;
	private final TurboIDExperimentType expType;

	private TurboID_Channel_Ori(Replicate replicate, TurboIDExperimentType expType) {
		this.replicate = replicate;
		this.expType = expType;
	}

	public Replicate getReplicate() {
		return replicate;
	}

	public static TurboID_Channel_Ori getByName(String name) {
		for (final TurboID_Channel_Ori value : values()) {
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
