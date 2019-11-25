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
	private final boolean useDistributedIntensity = true;
	private final boolean useNormalizedIntensity = false;

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

	public void run(boolean onlyTM, boolean onlyNonTM) throws IOException {
		for (final TurboIDExperimentType bait : TurboIDExperimentType.getBaits()) {
			String infix = "norm";
			if (useDistributedIntensity) {
				infix = "distr";
			}
			baitFile = new File(folder.getAbsoluteFile() + File.separator + bait.name() + "_" + infix + "_" + BAIT);
			preyFile = new File(folder.getAbsoluteFile() + File.separator + bait.name() + "_" + infix + "_" + PREY);
			interactionFile = new File(
					folder.getAbsoluteFile() + File.separator + bait.name() + "_" + infix + "_" + INTERACTION);
			writeBaitFile(bait);
			writePreyFile(onlyTM, onlyNonTM, bait);
			writeInteractionFile(onlyTM, onlyNonTM, bait);

		}
	}

	private void writeBaitFile(TurboIDExperimentType bait) throws IOException {
		FileWriter baitFW = null;

		try {
			baitFW = new FileWriter(baitFile);

			// applyLogs();
			this.controlExperiment.dealWithInfinities();
			this.testExperiment.dealWithInfinities();

			// control
			String testOrControl = "C";
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
//			}
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

	private void writePreyFile(boolean onlyTM, boolean onlyNonTM, TurboIDExperimentType bait) throws IOException {
		FileWriter preyFW = null;

		try {
			preyFW = new FileWriter(preyFile);
			final Set<String> accs = new THashSet<String>();
			accs.addAll(controlExperiment.keySet());
			accs.addAll(testExperiment.keySet());
			final UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(
					TurboIDDataAnalysis.uniprotReleasesFolder, true);
			final Map<String, Entry> annotatedProteins = uplr.getAnnotatedProteins(null, accs);
			final Set<String> preys = new THashSet<String>();
			for (final ProteinFromTurboID protein : controlExperiment.values()) {
				if (onlyTM && !protein.isTransmembrane()) {
					continue;
				}
				if (onlyNonTM && protein.isTransmembrane()) {
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
			for (final ProteinFromTurboID protein : testExperiment.values()) {
				if (onlyTM && !protein.isTransmembrane()) {
					continue;
				}
				if (onlyNonTM && protein.isTransmembrane()) {
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
			}
		}
	}

	private void writeInteractionFile(boolean onlyTM, boolean onlyNonTM, TurboIDExperimentType bait)
			throws IOException {
		FileWriter interactionFW = null;

		try {
			interactionFW = new FileWriter(interactionFile);
			// control
//			for (final TurboIDExperimentType bait : TurboIDExperimentType.getBaits()) {
			final String baitName = getBaitName(bait.name());
//				if (baitName == null) {
//					continue;
//				}
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
							final String preyName = protein.getGene();
							if (preyName.equals("Tbca") && baitName.equals("Sun1")) {
								log.info("asdf");
							}
							final int spc = protein.getSpc(replicate);
							double intensity = 0.0;
							if (useDistributedIntensity) {
								intensity = protein.getDistributedIntensitiesWithNormalizedIntensities().get(channel)
										* FACTOR;
							}
							if (useNormalizedIntensity) {
								intensity = protein.getNormalizedIntensities(bait).get(channel) * FACTOR;
							}
							if (!Double.isNaN(intensity) && Double.compare(0.0, intensity) != 0) {
								interactionFW
										.write(ipName + "\t" + baitName + "\t" + preyName + "\t" + intensity + "\n");
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
					if (channel.getExpType() != bait) {
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
							final String preyName = protein.getGene();
							final int spc = protein.getSpc(replicate);
							double intensity = 0.0;
							if (useDistributedIntensity) {
								intensity = protein.getDistributedIntensitiesWithNormalizedIntensities().get(channel)
										* FACTOR;
							}
							if (useNormalizedIntensity) {
								intensity = protein.getNormalizedIntensities(bait).get(channel) * FACTOR;
							}
							if (!Double.isNaN(intensity) && Double.compare(0.0, intensity) != 0) {
								interactionFW
										.write(ipName + "\t" + baitName + "\t" + preyName + "\t" + intensity + "\n");
							}
						}
					}
				}
			}
//			}

		} finally {
			if (interactionFW != null) {
				interactionFW.close();
			}
		}
	}

}
