package edu.scripps.yates.nucleome.turboID.panther;

import java.io.File;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import javax.xml.bind.JAXBException;

import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;

import edu.scripps.yates.annotations.panther.overrepresentation.PantherXMLOverrepresentationResultParser;
import edu.scripps.yates.annotations.panther.overrepresentation.PantherXMLOverrerepresentationResultUtil;
import edu.scripps.yates.annotations.panther.overrepresentation.xml.Overrepresentation;
import edu.scripps.yates.annotations.panther.overrepresentation.xml.Overrepresentation.Group;
import edu.scripps.yates.annotations.panther.overrepresentation.xml.Overrepresentation.Group.Result;
import edu.scripps.yates.annotations.panther.overrepresentation.xml.Overrepresentation.Group.Result.InputList;
import edu.scripps.yates.annotations.panther.overrepresentation.xml.Overrepresentation.Group.Result.Term;
import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.utilities.annotations.uniprot.UniprotEntryUtil;
import edu.scripps.yates.utilities.annotations.uniprot.xml.Entry;
import edu.scripps.yates.utilities.files.FileUtils;
import edu.scripps.yates.utilities.strings.StringUtils;
import edu.scripps.yates.utilities.util.Pair;

public class PantherResultsAnalyzer {
	private final static Logger log = Logger.getLogger(PantherResultsAnalyzer.class);
	private static final File resultsFolderProteomeAsReference = new File(
			"C:\\Users\\salvador\\Desktop\\4D_Nucleome\\TurboID\\input\\optimized_params\\slimGO\\fisher_FDR_genome_ref");
	private static final File resultsFolderQuantifiedAsReference = new File(
			"C:\\Users\\salvador\\Desktop\\4D_Nucleome\\TurboID\\input\\optimized_params\\slimGO\\fisher_FDR_quant_ref");
	private final UniprotProteinLocalRetriever uplr;
	private final static File uniprotFolder = new File("C:\\Users\\salvador\\Desktop\\uniprotKB");

	public PantherResultsAnalyzer() {
		this.uplr = new UniprotProteinLocalRetriever(uniprotFolder, true);
	}

	public static void main(String[] args) {
		PantherResultsAnalyzer p = new PantherResultsAnalyzer();
		try {
			p.run(resultsFolderProteomeAsReference, "Whole mouse genome");
		} catch (final JAXBException e) {
			e.printStackTrace();
		} catch (final IOException e) {
			e.printStackTrace();
		}

		p = new PantherResultsAnalyzer();
		try {
			p.run(resultsFolderQuantifiedAsReference, "2574 proteins quantified");
		} catch (final JAXBException e) {
			e.printStackTrace();
		} catch (final IOException e) {
			e.printStackTrace();
		}
	}

	private void run(File resultsFolder, String reference) throws JAXBException, IOException {
		final File excelFile = new File(resultsFolder.getAbsolutePath() + File.separator + "Panther_results.xlsx");
		if (excelFile.exists()) {
			excelFile.delete();
		}
		final List<File> resultsFiles = getResultsFiles(resultsFolder);
		final List<Overrepresentation> pantherResults = parsePantherResults(resultsFiles);
		for (final Overrepresentation overrepresentation : pantherResults) {
			final File file = printTableFile(resultsFolder, overrepresentation, reference);
			final String sheetName = FilenameUtils.getBaseName(file.getAbsolutePath());
			FileUtils.separatedValuesToXLSX(file.getAbsolutePath(), excelFile.getAbsolutePath(), "\t", sheetName);
			log.info("File " + file.getAbsolutePath() + " added to excel in sheet  " + sheetName);

			final File fileForR = printRFile(new File(resultsFolder.getAbsolutePath() + File.separator + "R"),
					sheetName, overrepresentation);
			log.info("File " + fileForR.getAbsolutePath() + " for R");
			file.delete();
		}
		log.info("Excel file at: " + excelFile.getAbsolutePath());
	}

	private File printRFile(File folder, String name, Overrepresentation overrepresentation) throws IOException {
		if (!folder.exists()) {
			folder.mkdirs();
		}
		final FileWriter fw = new FileWriter(folder.getAbsolutePath() + File.separator + name + "_edges.txt");
		fw.write("from\tto\n");
		final FileWriter fw2 = new FileWriter(folder.getAbsolutePath() + File.separator + name + "_vertices.txt");
		fw2.write("name\tfdr\n");
		fw2.write("root\t0\n");
		final List<Group> groups = overrepresentation.getGroup();
		for (final Group group : groups) {
			final List<Result> results = PantherXMLOverrerepresentationResultUtil.getOverrepresentedResults(group);
			if (results.isEmpty()) {
				continue;
			}
			String parent = null;
			String children = null;
			int previousLevel = -1;
			if (results.size() > 1) {
				for (int i = 0; i < results.size(); i++) {
					final Result result = results.get(i);
					final int level = result.getTerm().getLevel();
					final String label = result.getTerm().getLabel();
					// in fw2 we print the label and the fdr
					fw2.write(label + "\t" + -Math.log10(result.getInputList().getFdr()) + "\n");
					if (previousLevel == -1) {
						fw.write("root\t" + label + "\n");
						parent = label;
					} else {
						if (previousLevel < level) {
							children = label;
							fw.write(parent + "\t" + children + "\n");
						} else {
							parent = null;
							// previous level was the same or higher, so now we need to find the parent
							for (int j = i - 1; j >= 0; j--) {
								final int level2 = results.get(j).getTerm().getLevel();
								if (level2 < level) {
									parent = results.get(j).getTerm().getLabel();
									children = label;
									break;
								}
							}
							if (parent != null) {
								fw.write(parent + "\t" + children + "\n");
							} else {
								parent = label;
							}
						}
					}
					previousLevel = level;
					if (children != null) {
						parent = children;
					}
				}
			} else {
				fw.write("root\t" + results.get(0).getTerm().getLabel() + "\n");
				// in fw2 we print the label and the fdr
				fw2.write(results.get(0).getTerm().getLabel() + "\t"
						+ -Math.log10(results.get(0).getInputList().getFdr()) + "\n");
			}

		}

		fw.close();
		fw2.close();
		return folder;
	}

	private File printTableFile(File resultsFolder, Overrepresentation overrepresentation, String reference2)
			throws IOException {
		final File file = new File(
				resultsFolder.getAbsolutePath() + File.separator + getFileNameFromResults(overrepresentation) + ".txt");
		final FileWriter fw = new FileWriter(file);

		// header
		fw.write("Input file:\t" + PantherXMLOverrerepresentationResultUtil.getInputList(overrepresentation) + "\n");
		fw.write("Bait:\t" + getBaitName(overrepresentation) + "\n");

		fw.write("Reference:\t" + reference2 + "\n");
		fw.write("Organism selected for input list:\t"
				+ PantherXMLOverrerepresentationResultUtil.getInputFileOrganism(overrepresentation) + "\n");
		fw.write("Test type:\t" + overrepresentation.getTestType() + "\n");
		fw.write("p-value Correction:\t" + overrepresentation.getCorrection() + "\n");
		fw.write("Annotation type:\t" + overrepresentation.getAnnotationType() + "\n");
		fw.write("Panther version:\t" + overrepresentation.getDataVersionReleaseDate() + "\n");
		fw.write("Terms overrepresented:\t"
				+ PantherXMLOverrerepresentationResultUtil.getNumTermsOverrepresented(overrepresentation) + "\n");
		fw.write("Terms underrepresented:\t"
				+ PantherXMLOverrerepresentationResultUtil.getNumTermsUnderrepresented(overrepresentation) + "\n");
		fw.write("Groups overrepresented:\t"
				+ PantherXMLOverrerepresentationResultUtil.getNumGroupTermsOverrepresented(overrepresentation) + "\n");
		fw.write("Groups underrepresented:\t"
				+ PantherXMLOverrerepresentationResultUtil.getNumGroupTermsUnderrepresented(overrepresentation) + "\n");
		fw.write("\n\n");
		fw.write(
				"level\tGO\tlabel\tNumber in reference\tNumber in list\tExpected\tPlus_Minus\tp-value\tfold_enrichment\tFDR\tProteins\tGenes\n");

		for (final Group group : overrepresentation.getGroup()) {
			final List<Result> results = group.getResult();
			for (final Result result : results) {
				final Term term = result.getTerm();
				fw.write(term.getLevel() + "\t");
				fw.write(term.getId() + "\t");
				fw.write(term.getLabel() + "\t");
				fw.write(result.getNumberInReference() + "\t");
				final InputList inputList = result.getInputList();
				fw.write(inputList.getNumberInList() + "\t");
				fw.write(inputList.getExpected() + "\t");
				fw.write(inputList.getPlusMinus() + "\t");
				fw.write(inputList.getPValue() + "\t");
				fw.write(inputList.getFoldEnrichment() + "\t");
				fw.write(inputList.getFdr() + "\t");
				if (inputList.getMappedIdList() != null) {
					final List<String> mappedIds = inputList.getMappedIdList().getMappedId();
					Collections.sort(mappedIds);
					final List<Object> list = new ArrayList<Object>();
					list.addAll(mappedIds);
					fw.write(StringUtils.getSeparatedValueStringFromChars(list, ",") + "\t");
					final Map<String, Entry> annotatedProteins = uplr.getAnnotatedProteins(null, mappedIds);

					final StringBuilder sb = new StringBuilder();
					for (final String acc : mappedIds) {
						if (annotatedProteins.containsKey(acc)) {
							final List<Pair<String, String>> geneNames = UniprotEntryUtil
									.getGeneName(annotatedProteins.get(acc), true, true);
							if (!geneNames.isEmpty()) {

								final String gene = geneNames.get(0).getFirstelement();
								sb.append(gene);
							} else {
								log.info("No gene for protein " + acc);
								sb.append("N/A");
							}
						}
						sb.append(",");
					}
					fw.write(sb.toString());
				}
				fw.write("\n");
			}
		}
		fw.close();
		log.info("File written at: " + file.getAbsolutePath());
		return file;
	}

	private String getFileNameFromResults(Overrepresentation overrepresentation) {
		final String baitName = getBaitName(overrepresentation);
		final String[] split = overrepresentation.getAnnotationType().split(" ");
		final String goType = split[2].charAt(0) + "" + split[3].charAt(0);
		return baitName + "_" + goType;
	}

	private String getBaitName(Overrepresentation overrepresentation) {
		return FilenameUtils.getBaseName(PantherXMLOverrerepresentationResultUtil.getInputList(overrepresentation))
				.split("_")[0];
	}

	private List<Overrepresentation> parsePantherResults(List<File> resultsFiles) throws JAXBException {
		final List<Overrepresentation> ret = new ArrayList<Overrepresentation>();
		for (final File pantherFile : resultsFiles) {
			final PantherXMLOverrepresentationResultParser parser = new PantherXMLOverrepresentationResultParser(
					pantherFile);
			ret.add(parser.getOverrepresentation());
		}
		return ret;
	}

	private List<File> getResultsFiles(File resultsFolder) {
		final File[] array = resultsFolder.listFiles(new FilenameFilter() {

			@Override
			public boolean accept(File dir, String name) {
				if (name.contains("analysis.xml")) {
					return true;
				}
				return false;
			}
		});
		return Arrays.asList(array);
	}
}
