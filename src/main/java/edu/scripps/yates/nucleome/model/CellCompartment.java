package edu.scripps.yates.nucleome.model;

/**
 * NE, N, C
 * 
 * @author Salva
 *
 */
public enum CellCompartment {
	NE, CML, N, CMH_NEW, CMH_OLD;

	public static CellCompartment[] getCellCompartmentsForXiDatasets() {
		final CellCompartment[] ret = new CellCompartment[4];
		ret[0] = CML;
		ret[1] = CMH_NEW;
		ret[2] = CMH_OLD;
		ret[3] = NE;
		return ret;
	}
}
