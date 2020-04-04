package edu.scripps.yates.nucleome.turboID;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.utilities.annotations.uniprot.UniprotEntryUtil;
import edu.scripps.yates.utilities.annotations.uniprot.xml.Entry;
import gnu.trove.set.hash.THashSet;

/**
 * This class has method to generate the required files to use in as input with
 * SAINTExpress
 * 
 * @author salvador
 *
 */
public class SAINTExpress {
	private final static String BAIT = "bait.txt";
	private final static String PREY = "prey.txt";
	private final static String INTERACTION = "interaction.txt";
	private final File folder;
	private final TurboIDExperiment controlExperiment;
	private final TurboIDExperiment testExperiment;
	private File baitFile;
	private File preyFile;
	private File interactionFile;

	private final static Logger log = Logger.getLogger(SAINTExpress.class);
	private static final double FACTOR = 1000000;
	//////////////////////////
	//
	private boolean useDistributedIntensity = true;
	private boolean useNormalizedIntensity = false;
	private boolean runAllBaitsTogether = false;

	//
	///////////////////////
	public SAINTExpress(File folder, TurboIDExperiment controlExperiment, TurboIDExperiment testExperiment)
			throws IOException {
		this.folder = folder;
		this.controlExperiment = controlExperiment;
		this.testExperiment = testExperiment;
		if (!folder.exists()) {
			folder.mkdirs();
		}

	}

	public void run(boolean onlyTM, boolean onlyNonTM, boolean onlyComplete, int minReplicates) throws IOException {
		if (runAllBaitsTogether) {
			runBaitsTogether(onlyTM, onlyNonTM, onlyComplete, minReplicates);

		} else {
			runBaitsIndependently(onlyTM, onlyNonTM, onlyComplete, minReplicates);

		}
	}

	public void runBaitsTogether(boolean onlyTM, boolean onlyNonTM, boolean onlyComplete, int minReplicates)
			throws IOException {
		String infix = "norm";
		if (useDistributedIntensity) {
			infix = "distr";
		}
		String suffix = "";
		if (controlExperiment == null) {
			suffix = "_TIDonlyCTRL";
		}
		baitFile = new File(folder.getAbsoluteFile() + File.separator + "All_Baits_" + infix + suffix + "_" + BAIT);
		preyFile = new File(folder.getAbsoluteFile() + File.separator + "All_Baits_" + infix + suffix + "_" + PREY);
		interactionFile = new File(
				folder.getAbsoluteFile() + File.separator + "All_Baits_" + infix + suffix + "_" + INTERACTION);
		boolean append = false;
		writePreyFile(onlyTM, onlyNonTM, onlyComplete, minReplicates, append);
		int i = 0;
		for (final TurboIDExperimentType bait : TurboIDExperimentType.getBaits()) {
			append = false;
			if (i > 0) {
				append = true;
			}
			boolean writeControl = true;
			// if we use the option to use all baits together, and controlExperiment is
			// null, that is, the control is the TBID channel, the control baits lines will
			// be repeated after the first bait, because all baits share the same control

			if (controlExperiment == null && i > 0) {
				writeControl = false;
			}
			writeBaitFile(writeControl, bait, append);
			writeInteractionFile(onlyTM, onlyNonTM, onlyComplete, minReplicates, bait, append);
			i++;
		}
	}

	public void runBaitsIndependently(boolean onlyTM, boolean onlyNonTM, boolean onlyComplete, int minReplicates)
			throws IOException {
		final boolean append = false;
		String infix = "norm";
		if (useDistributedIntensity) {
			infix = "distr";
		}
		String suffix = "";
		if (controlExperiment == null) {
			suffix = "_TIDonlyCTRL";
		}
		for (final TurboIDExperimentType bait : TurboIDExperimentType.getBaits()) {
			baitFile = new File(
					folder.getAbsoluteFile() + File.separator + bait.name() + "_" + infix + suffix + "_" + BAIT);
			preyFile = new File(
					folder.getAbsoluteFile() + File.separator + bait.name() + "_" + infix + suffix + "_" + PREY);
			interactionFile = new File(
					folder.getAbsoluteFile() + File.separator + bait.name() + "_" + infix + suffix + "_" + INTERACTION);
			writeBaitFile(true, bait, append);
			writePreyFile(onlyTM, onlyNonTM, onlyComplete, minReplicates, append);
			writeInteractionFile(onlyTM, onlyNonTM, onlyComplete, minReplicates, bait, append);

		}
	}

	private void writeBaitFile(boolean writeControl, TurboIDExperimentType bait, boolean append) throws IOException {
		FileWriter baitFW = null;

		try {
			baitFW = new FileWriter(baitFile, append);

			// applyLogs();
			String testOrControl = "C";
			this.testExperiment.dealWithInfinities();
			if (writeControl) {

				if (controlExperiment != null) {
					this.controlExperiment.dealWithInfinities();

					// control

//			for (final TurboIDExperimentType bait : TurboIDExperimentType.getBaits()) {
					for (final Replicate replicate : Replicate.values(controlExperiment.getFraction())) {
						// do not use Arerun
						if (replicate == Replicate.Arerun) {
							continue;
						}
						for (final TurboID_Channel_Norm channel : TurboID_Channel_Norm.values()) {
							if (channel.getExpType() != bait) {
								continue;
							}
							if (channel.getReplicate() == replicate) {
								final String ipName = replicate.name() + "_" + channel.name();
								final String baitName = getBaitName(bait.name()) + "_C";
								if (baitName != null) {
									baitFW.write(ipName + "\t" + baitName + "\t" + testOrControl + "\n");
								}
							}
						}
					}
				} else {
					// if control is null, then use the TURBO_ID_ONLY_AS_CONTROL from the test
					// experiment as control
					for (final Replicate replicate : Replicate.values(testExperiment.getFraction())) {
						// do not use Arerun
						if (replicate == Replicate.Arerun) {
							continue;
						}
						for (final TurboID_Channel_Norm channel : TurboID_Channel_Norm.values()) {
							if (channel.getExpType() != TurboIDExperimentType.TURBOID_ONLY) {
								continue;
							}
							if (channel.getReplicate() == replicate) {
								final String ipName = replicate.name() + "_" + channel.name();
								final String baitName = getBaitName(bait.name()) + "_C";
								if (baitName != null) {
									baitFW.write(ipName + "\t" + baitName + "\t" + testOrControl + "\n");
								}
							}
						}

					}
				}
			}
			// test
			testOrControl = "T";
//			for (final TurboIDExperimentType bait : TurboIDExperimentType.getBaits()) {
			for (final Replicate replicate : Replicate.values(testExperiment.getFraction())) {
				// do not use Arerun
				if (replicate == Replicate.Arerun) {
					continue;
				}
				for (final TurboID_Channel_Norm channel : TurboID_Channel_Norm.values()) {
					if (channel.getExpType() != bait) {
						continue;
					}
					// if control is null, then use the TURBO_ID_ONLY_AS_CONTROL, so here, skip them
					if (controlExperiment == null && channel.getExpType() == TurboIDExperimentType.TURBOID_ONLY) {
						continue;
					}
					if (channel.getReplicate() == replicate) {
						final String ipName = replicate.name() + "_" + channel.name();

						final String baitName = getBaitName(bait.name()) + "_T";
						;
						if (baitName != null) {
							baitFW.write(ipName + "\t" + baitName + "\t" + testOrControl + "\n");
						}
					}
				}
			}
//			}

		} finally {
			if (baitFW != null) {
				baitFW.close();
				log.info("Bait file written at: " + interactionFile.getAbsolutePath());
			}
		}
	}

	private String getBaitName(String name) {
		if (name.equals(TurboIDExperimentType.EMD.name())) {
			return "Emd";
		} else if (name.equals(TurboIDExperimentType.LBR.name())) {
			return "Lbr";
		} else if (name.equals(TurboIDExperimentType.MAN1.name())) {
			return "Man1";
		} else if (name.equals(TurboIDExperimentType.SUN1.name())) {
			return "Sun1";
		}
		return null;
	}

	private void writePreyFile(boolean onlyTM, boolean onlyNonTM, boolean onlyComplete, int minReps, boolean append)
			throws IOException {
		FileWriter preyFW = null;

		try {
			preyFW = new FileWriter(preyFile, append);
			final Set<String> accs = new THashSet<String>();
			if (controlExperiment != null) {
				accs.addAll(controlExperiment.keySet());
			}
			accs.addAll(testExperiment.keySet());
			final UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(
					TurboIDDataAnalysis.uniprotReleasesFolder, true);
			final Map<String, Entry> annotatedProteins = uplr.getAnnotatedProteins(null, accs);
			final Set<String> preys = new THashSet<String>();
			if (controlExperiment != null) {
				for (final ProteinFromTurboID protein : controlExperiment.values()) {
					if (onlyTM && !protein.isTransmembrane()) {
						continue;
					}
					if (onlyNonTM && protein.isTransmembrane()) {
						continue;
					}
					if (onlyComplete && !protein.isComplete(minReps, controlExperiment.getFraction())) {
						continue;
					}
					final String preyName = protein.getGene();
					if ("N/A".equals(preyName)) {
						log.info(protein.getAcc() + " ignoring it");
						continue;
					}
					if (annotatedProteins.containsKey(protein.getAcc())) {
						if (preys.contains(preyName)) {
							continue;
						}
						preys.add(preyName);
						preyFW.write(preyName + "\t");
						final Entry entry = annotatedProteins.get(protein.getAcc());
						final int length = UniprotEntryUtil.getProteinSequence(entry).length();
						preyFW.write(length + "\t");
						preyFW.write(preyName + "\n");
					} else {
						log.debug("Protein " + preyName + " ignored because it is not found in Uniprot");
					}

				}
			}
			for (final ProteinFromTurboID protein : testExperiment.values()) {
				if (onlyTM && !protein.isTransmembrane()) {
					continue;
				}
				if (onlyNonTM && protein.isTransmembrane()) {
					continue;
				}
				if (onlyComplete && !protein.isComplete(minReps, testExperiment.getFraction())) {
					continue;
				}
				if (annotatedProteins.containsKey(protein.getAcc())) {
					if (preys.contains(protein.getGene())) {
						continue;
					}
					preys.add(protein.getGene());
					preyFW.write(protein.getGene() + "\t");
					final Entry entry = annotatedProteins.get(protein.getAcc());
					final int length = UniprotEntryUtil.getProteinSequence(entry).length();
					preyFW.write(length + "\t");
					preyFW.write(protein.getGene() + "\n");
				} else {
					log.debug("Protein " + protein.getGene() + " ignored because it is not found in Uniprot");
				}

			}
		} finally {
			if (preyFW != null) {
				preyFW.close();
				log.info("Prey file written at: " + preyFile.getAbsolutePath());
			}
		}
	}

	private void writeInteractionFile(boolean onlyTM, boolean onlyNonTM, boolean onlyComplete, int minReps,
			TurboIDExperimentType bait, boolean append) throws IOException {
		FileWriter interactionFW = null;

		try {
			interactionFW = new FileWriter(interactionFile, append);
			// control
//			for (final TurboIDExperimentType bait : TurboIDExperimentType.getBaits()) {
			final String baitName = getBaitName(bait.name());
//				if (baitName == null) {
//					continue;
//				}
			if (controlExperiment != null) {
				for (final Replicate replicate : Replicate.values(controlExperiment.getFraction())) {
					// do not use Arerun
					if (replicate == Replicate.Arerun) {
						continue;
					}
					for (final TurboID_Channel_Norm channel : TurboID_Channel_Norm.values()) {
						if (channel.getExpType() != bait) {
							continue;
						}

						final String ipName = replicate.name() + "_" + channel.name();
						if (channel.getReplicate() == replicate) {
							for (final ProteinFromTurboID protein : controlExperiment.values()) {
								if (onlyTM && !protein.isTransmembrane()) {
									continue;
								}
								if (onlyNonTM && protein.isTransmembrane()) {
									continue;
								}
								if (onlyComplete && !protein.isComplete(minReps, controlExperiment.getFraction())) {
									continue;
								}
								final String preyName = protein.getGene();
								if ("N/A".equals(preyName)) {
									log.info(protein.getAcc() + " ignoring it");
									continue;
								}

								final int spc = protein.getSpc(replicate);
								double intensity = 0.0;
								if (useDistributedIntensity) {
									intensity = protein.getDistributedIntensitiesWithNormalizedIntensities()
											.get(channel) * FACTOR;
								}
								if (useNormalizedIntensity) {
									intensity = protein.getNormalizedIntensities(bait).get(channel) * FACTOR;
								}
								if (intensity == 56689.85813151278) {
									log.info("asdf");
								}

								if (!Double.isNaN(intensity) && Double.compare(0.0, intensity) != 0) {
									interactionFW.write(
											ipName + "\t" + baitName + "\t" + preyName + "\t" + intensity + "\n");
									if (preyName.equals("Lmna") && baitName.equals("Emd")) {
										log.info(ipName + "\t" + baitName + "\t" + preyName + "\t" + intensity);
									}
								}
							}
						}
					}
				}
			}
//			}
			// test
//			for (final TurboIDExperimentType bait : TurboIDExperimentType.getBaits()) {
//				  baitName = getBaitName(bait.name());
//				if (baitName == null) {
//					continue;
//				}
			for (final Replicate replicate : Replicate.values(testExperiment.getFraction())) {
				// do not use Arerun
				if (replicate == Replicate.Arerun) {
					continue;
				}
				for (final TurboID_Channel_Norm channel : TurboID_Channel_Norm.values()) {
					if (controlExperiment != null && channel.getExpType() != bait) {
						continue;
					}
					// if there is no control experiment, here we print only the ones that are like
					// the bait OR the ones that are from TURBO_ID_ONLY
					if (controlExperiment == null && channel.getExpType() != bait
							&& channel.getExpType() != TurboIDExperimentType.TURBOID_ONLY) {
						continue;
					}
					final String ipName = replicate.name() + "_" + channel.name();
					if (channel.getReplicate() == replicate) {
						for (final ProteinFromTurboID protein : testExperiment.values()) {
							if (onlyTM && !protein.isTransmembrane()) {
								continue;
							}
							if (onlyNonTM && protein.isTransmembrane()) {
								continue;
							}
							if (onlyComplete && !protein.isComplete(minReps, testExperiment.getFraction())) {
								continue;
							}
							final String preyName = protein.getGene();
							if ("N/A".equals(preyName)) {
								log.info(protein.getAcc() + " ignoring it");
								continue;
							}
							final int spc = protein.getSpc(replicate);
							double intensity = 0.0;
							if (useDistributedIntensity) {
								intensity = protein.getDistributedIntensitiesWithNormalizedIntensities().get(channel)
										* FACTOR;
							}

							if (useNormalizedIntensity) {
								intensity = protein.getNormalizedIntensities(channel.getExpType()).get(channel)
										* FACTOR;
							}
							if (!Double.isNaN(intensity) && Double.compare(0.0, intensity) != 0) {
								interactionFW
										.write(ipName + "\t" + baitName + "\t" + preyName + "\t" + intensity + "\n");
								if (preyName.equals("Lmna") && baitName.equals("Emd")) {
									log.info(ipName + "\t" + baitName + "\t" + preyName + "\t" + intensity);
								}
							}
						}
					}
				}
			}
//			}

		} finally {
			if (interactionFW != null) {
				interactionFW.close();
				log.info("Interaction file written at: " + interactionFile.getAbsolutePath());
			}
		}
	}

	public boolean isUseDistributedIntensity() {
		return useDistributedIntensity;
	}

	public void setUseDistributedIntensity(boolean useDistributedIntensity) {
		this.useDistributedIntensity = useDistributedIntensity;
	}

	public boolean isUseNormalizedIntensity() {
		return useNormalizedIntensity;
	}

	public void setUseNormalizedIntensity(boolean useNormalizedIntensity) {
		this.useNormalizedIntensity = useNormalizedIntensity;
	}

	public void setRunAllBaitsTogether(boolean runAllBaitsTogether) {
		this.runAllBaitsTogether = runAllBaitsTogether;
	}
}
