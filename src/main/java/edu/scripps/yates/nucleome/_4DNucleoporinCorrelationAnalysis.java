package edu.scripps.yates.nucleome;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.poi.hssf.util.CellReference;

import edu.scripps.yates.excel.ExcelReader;
import edu.scripps.yates.nucleome.model.CellType;

public class _4DNucleoporinCorrelationAnalysis {
	private final static String TAB = "\t";

	public static void main(String[] args) {
		File inputFile = new File(
				"Z:\\share\\Salva\\data\\4D_Nucleome\\SwissProt_1FDR\\output\\Nucleoporin_correlations.xlsx");
		ExcelReader reader;
		try {
			reader = new ExcelReader(inputFile, 0, 0);
			List<String> proteinList = new ArrayList<String>();
			Map<String, Double> correlations = new HashMap<String, Double>();
			final int lastRow = 54;
			for (int row1 = 0; row1 < lastRow; row1++) {
				String gene1 = reader.getStringValue(0, row1, "C");
				if (!gene1.equals("Gene")) {
					if (!proteinList.contains(gene1)) {
						proteinList.add(gene1);
					}
					String cellTypeString = reader.getStringValue(0, row1, "A");
					if (cellTypeString != null) {
						CellType cellType = CellType.valueOf(cellTypeString);

						String colBEGIN1 = getColBeginNSAF(cellType);
						String colEND1 = getColEndNSAF(cellType);
						double[] array1 = getArray(reader, row1, colBEGIN1, colEND1);
						for (int row2 = 0; row2 < lastRow; row2++) {
							// if (row2 == row1) {
							// continue;
							// }
							String gene2 = reader.getStringValue(0, row2, "C");
							if (!gene2.equals("Gene")) {

								String cellTypeString2 = reader.getStringValue(0, row2, "A");
								if (cellTypeString2 != null) {
									CellType cellType2 = CellType.valueOf(cellTypeString2);
									if (cellType != cellType2) {
										continue;
									}
									String comparisonKey = cellType + "-" + gene1 + "-" + gene2;

									String colBEGIN2 = getColBeginNSAF(cellType2);
									String colEND2 = getColEndNSAF(cellType2);

									double[] array2 = getArray(reader, row2, colBEGIN2, colEND2);
									double correlation = new PearsonsCorrelation().correlation(array1, array2);
									System.out.println(cellType + TAB + gene1 + TAB + gene2 + TAB + correlation);
									correlations.put(comparisonKey, correlation);
								}
							}
						}

					}
				}
			}
			for (CellType cellType : CellType.values()) {
				System.out.println(cellType);
				// write header
				System.out.print(TAB);
				for (String protein : proteinList) {
					System.out.print(protein + TAB);
				}
				System.out.println();
				for (int i = 0; i < proteinList.size(); i++) {
					String gene1 = proteinList.get(i);
					System.out.print(gene1 + TAB);
					for (int j = 0; j < proteinList.size(); j++) {
						String gene2 = proteinList.get(j);
						String comparisonKey = cellType + "-" + gene1 + "-" + gene2;
						Double corr = correlations.get(comparisonKey);
						System.out.print(corr + TAB);
					}
					System.out.println();
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private static String getColEndNSAF(CellType cellType) {
		if (cellType == CellType.A) {
			return "N";
		} else if (cellType == CellType.U) {
			return "R";
		} else if (cellType == CellType.M) {
			return "P";
		}
		throw new IllegalArgumentException(cellType + " is a not supported cellType");
	}

	private static String getColEndSPC(CellType cellType) {
		if (cellType == CellType.A) {
			return "Q";
		} else if (cellType == CellType.U) {
			return "M";
		} else if (cellType == CellType.M) {
			return "K";
		}
		throw new IllegalArgumentException(cellType + " is a not supported cellType");
	}

	private static String getColBeginNSAF(CellType cellType) {
		if (cellType == CellType.A) {
			return "L";
		} else if (cellType == CellType.U) {
			return "N";
		} else if (cellType == CellType.M) {
			return "M";
		}
		throw new IllegalArgumentException(cellType + " is a not supported cellType");
	}

	private static String getColBeginSPC(CellType cellType) {
		return "I";
	}

	private static double[] getArray(ExcelReader reader, int row, String colBEGIN, String colEND) {
		int start = CellReference.convertColStringToIndex(colBEGIN);
		int end = CellReference.convertColStringToIndex(colEND);
		double[] ret = new double[end - start + 1];
		int i = 0;
		for (int col = start; col <= end; col++) {
			ret[i++] = Double.valueOf(reader.getNumberValue(0, row, col));
		}
		return ret;
	}

}
