package edu.scripps.yates.nucleome.turboID;

import org.apache.log4j.Logger;

import gnu.trove.map.hash.THashMap;

public class IPExperiment extends THashMap<String, ProteinFromIP> {
	private final static Logger log = Logger.getLogger(IPExperiment.class);

//	public File exportToFileForClustering(File outputFile, boolean onlyTM, boolean onlyNonTM,
//			boolean ignoreMissingValues, boolean useNormalizedIntensities) {
//
////		applyLogs();
//		dealWithInfinities();
//
//		FileWriter fw = null;
//		try {
//			fw = new FileWriter(outputFile);
//			// header
//			fw.write("ACC\t");
//			fw.write("Gene\t");
//			fw.write("Transmembrane\t");
//			// spc
//			for (final Replicate replicate : Replicate.values()) {
//				fw.write("SPC " + replicate.name() + "\t");
//			}
//			StringBuilder sb = new StringBuilder();
//			for (final TurboIDExperimentType turboIDExperimentType : TurboIDExperimentType.values()) {
//				if (turboIDExperimentType != TurboIDExperimentType.MIX) {
//					if (!"".equals(sb.toString())) {
//						sb.append("\t");
//					}
//					sb.append(turboIDExperimentType);
//					if (turboIDExperimentType == TurboIDExperimentType.TURBOID_ONLY) {
//
//						for (final Replicate replicate : Replicate.values()) {
//							for (final Channel_Norm channel : Channel_Norm.values()) {
//								if (channel.getReplicate() == replicate) {
//									if (channel.name().contains("10") || channel.name().contains("11")) {
//
//										if (!"".equals(sb.toString())) {
//											sb.append("\t");
//										}
//										sb.append(channel.getReplicate() + "-" + channel.name());
//									}
//								}
//							}
//
//						}
//					}
//				}
//			}
//			fw.write(sb.toString());
//			fw.write("\n");
//			for (final ProteinFromIP protein : values()) {
//				if (onlyTM && !protein.isTransmembrane()) {
//					continue;
//				}
//				if (onlyNonTM && protein.isTransmembrane()) {
//					continue;
//				}
//				try {
//					fw.write(protein.getAcc() + "\t");
//					fw.write(protein.getGene() + "\t");
//					fw.write(protein.isTransmembrane() + "\t");
//					// spc
//					for (final Replicate replicate : Replicate.values()) {
//						fw.write(protein.getSpc(replicate) + "\t");
//					}
//					sb = new StringBuilder();
//					for (final TurboIDExperimentType turboIDExperimentType : TurboIDExperimentType.values()) {
//						if (turboIDExperimentType != TurboIDExperimentType.MIX) {
//							double avg = Double.NaN;
//							if (useNormalizedIntensities) {
//								avg = protein.getAverageFromNormalized(turboIDExperimentType, ignoreMissingValues);
//							} else {
//								avg = protein.getAverageFromOriginal(turboIDExperimentType, ignoreMissingValues);
//							}
//							if (!"".equals(sb.toString())) {
//								sb.append("\t");
//							}
//							sb.append(avg);
//							if (turboIDExperimentType == TurboIDExperimentType.TURBOID_ONLY) {
//								final TObjectDoubleHashMap<Channel_Norm> normalizedIntensities = protein
//										.getNormalizedIntensities(turboIDExperimentType);
//								for (final Replicate replicate : Replicate.values()) {
//									for (final Channel_Norm channel : Channel_Norm.values()) {
//										if (channel.getReplicate() == replicate) {
//											if (channel.name().contains("10") || channel.name().contains("11")) {
//												final double intensity = normalizedIntensities.get(channel);
//												if (!"".equals(sb.toString())) {
//													sb.append("\t");
//												}
//												sb.append(intensity);
//											}
//										}
//									}
//
//								}
//							}
//						}
//					}
//					fw.write(sb.toString());
//				} finally {
//					fw.write("\n");
//				}
//			}
//		} catch (final IOException e) {
//			e.printStackTrace();
//		} finally {
//			if (fw != null) {
//				try {
//					fw.close();
//				} catch (final IOException e) {
//				}
//			}
//		}
//		return outputFile;
//	}

//	private void applyLogs() {
//		for (final ProteinFromTurboID protein : values()) {
//			for (final TurboIDExperimentType expType : TurboIDExperimentType.values()) {
//				final TObjectDoubleHashMap<Channel_Norm> intensities = protein.getNormalizedIntensities(expType);
//				for (final Channel_Norm channel : intensities.keySet()) {
//					final double d = intensities.get(channel);
//					intensities.put(channel, Maths.log(d, 2));
//				}
//				final TObjectDoubleHashMap<Channel_Ori> intensities2 = protein.getOriginalIntensities(expType);
//				for (final Channel_Ori channel : intensities2.keySet()) {
//					final double d = intensities2.get(channel);
//					intensities2.put(channel, Maths.log(d, 2));
//				}
//			}
//		}
//	}

	private void dealWithInfinities() {
		for (final IPExperimentType expType : IPExperimentType.values()) {
			Double maxNorm = -Double.MAX_VALUE;
			Double minNorm = Double.MAX_VALUE;
			Double maxOri = -Double.MAX_VALUE;
			Double minOri = Double.MAX_VALUE;
			for (final ProteinFromIP protein : values()) {
				Double d = protein.getNormalizedIntensity(expType);

				if (Double.isFinite(d)) {
					if (d > maxNorm) {
						maxNorm = d;
					}
					if (d < minNorm) {
						minNorm = d;
					}
				}

				d = protein.getOriginalIntensity(expType);

				if (Double.isFinite(d)) {
					if (d > maxOri) {
						maxOri = d;
					}
					if (d < minOri) {
						minOri = d;
					}
				}

			}
			log.info("Maximum from normalized for " + expType + ": " + maxNorm);
			log.info("Minimum from normalized for " + expType + ": " + minNorm);
			log.info("Maximum from original for " + expType + ": " + maxOri);
			log.info("Minimum from original for " + expType + ": " + minOri);
			//
			for (final ProteinFromIP protein : values()) {
				final Double intensity = protein.getNormalizedIntensity(expType);

				if (Double.isInfinite(intensity)) {
					if (Double.compare(Double.NEGATIVE_INFINITY, intensity) == 0) {
						protein.addNormalizedIntensity(minNorm, expType);
					}
					if (Double.compare(Double.POSITIVE_INFINITY, intensity) == 0) {
						protein.addNormalizedIntensity(maxNorm, expType);
					}
				}

				final Double intensity2 = protein.getOriginalIntensity(expType);
				if (Double.isInfinite(intensity2)) {
					if (Double.compare(Double.NEGATIVE_INFINITY, intensity2) == 0) {
						protein.addOriginalIntensity(minOri, expType);
					}
					if (Double.compare(Double.POSITIVE_INFINITY, intensity2) == 0) {
						protein.addOriginalIntensity(maxOri, expType);
					}
				}

			}
		}
	}

}
