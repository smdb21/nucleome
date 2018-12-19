package edu.scripps.yates.nucleome.model;

import java.util.ArrayList;
import java.util.List;

/**
 * UREA or CARBONATE
 * 
 * 
 * @author Salva
 *
 */
public enum Wash {
	U1, U2, U3, U4,
	//
	//
	C1, C3,
	// C4,
	preC1, preU1, preU2, preU3, preU4,
	//

	// , WC
	;

	public static List<Wash> getWashesWithAllFractionsFromWashGroup(WashGroup washGroup) {
		final List<Wash> ret = new ArrayList<Wash>();
		switch (washGroup) {
		case U:
			ret.add(U1);
			ret.add(U2);
			ret.add(U3);
			ret.add(U4);

			break;
		case C:
			ret.add(C1);
			ret.add(C3);
			break;
		case Pre:
			ret.add(preC1);
			ret.add(preU1);
			ret.add(preU4);
			ret.add(preU2);
			ret.add(preU3);

			break;
		default:
			break;
		}
		return ret;
	}
}
