package edu.scripps.yates.nucleome;

import edu.scripps.yates.nucleome.model.Wash;
import edu.scripps.yates.nucleome.model.WashGroup;
import edu.scripps.yates.utilities.util.Pair;
import gnu.trove.map.hash.TObjectDoubleHashMap;

public class CalculatedScores {

	private final String proteinKey;
	private final TObjectDoubleHashMap<Wash> individualScoresWithNoCMHFraction = new TObjectDoubleHashMap<Wash>();
	private final TObjectDoubleHashMap<Wash> individualScoresWithCMHFraction = new TObjectDoubleHashMap<Wash>();
	private final TObjectDoubleHashMap<WashGroup> groupedScoresByWashGroup = new TObjectDoubleHashMap<WashGroup>();
	private final TObjectDoubleHashMap<Pair<Wash, Wash>> nuclearEnvelopeFoldChangesPreWashVSWash = new TObjectDoubleHashMap<Pair<Wash, Wash>>();
	private final TObjectDoubleHashMap<Pair<Wash, Wash>> nuclearEnvelopeFoldChangesPValuesPreWashVSWash = new TObjectDoubleHashMap<Pair<Wash, Wash>>();
	private Integer totalNuclearContentFractionSPCInUInSabyDataset;
	private Double nuclearContentEnrichmentScoreInUInSabyDataset;

	public CalculatedScores(String proteinKey) {
		this.proteinKey = proteinKey;
	}

	public void addIndividualScoreWithNoCMHFraction(double score, Wash wash) {
		individualScoresWithNoCMHFraction.put(wash, score);
	}

	public double getIndividualScoreWithNoCMHFraction(Wash wash) {
		if (!individualScoresWithNoCMHFraction.containsKey(wash)) {
			throw new IllegalArgumentException(wash + " not in the hash");
		}
		return individualScoresWithNoCMHFraction.get(wash);
	}

	public void addIndividualScoreWithCMHFraction(double score, Wash wash) {
		individualScoresWithCMHFraction.put(wash, score);
	}

	public double getIndividualScoreWithCMHFraction(Wash wash) {
		if (!individualScoresWithCMHFraction.containsKey(wash)) {
			throw new IllegalArgumentException(wash + " not in the hash");
		}
		return individualScoresWithCMHFraction.get(wash);
	}

	public void addGroupedScore(double score, WashGroup washGroup) {
		groupedScoresByWashGroup.put(washGroup, score);
	}

	public double getGroupedScore(WashGroup washGroup) {
		if (!groupedScoresByWashGroup.containsKey(washGroup)) {
			throw new IllegalArgumentException(washGroup + " not in the hash");
		}
		return groupedScoresByWashGroup.get(washGroup);
	}

	public void addNEFoldChange(double foldChange, Wash washPostWash, Wash washPreWash) {
		final Pair<Wash, Wash> pair = new Pair<Wash, Wash>(washPostWash, washPreWash);
		nuclearEnvelopeFoldChangesPreWashVSWash.put(pair, foldChange);
	}

	public double getNEFoldChange(Wash washPostWash, Wash washPreWash) {
		final Pair<Wash, Wash> pair = new Pair<Wash, Wash>(washPostWash, washPreWash);
		if (!nuclearEnvelopeFoldChangesPreWashVSWash.containsKey(pair)) {
			throw new IllegalArgumentException(pair + " not in the hash");
		}
		return nuclearEnvelopeFoldChangesPreWashVSWash.get(pair);
	}

	public void addNEFoldChangePValue(double pValue, Wash washPostWash, Wash washPreWash) {
		final Pair<Wash, Wash> pair = new Pair<Wash, Wash>(washPostWash, washPreWash);
		nuclearEnvelopeFoldChangesPValuesPreWashVSWash.put(pair, pValue);
	}

	public double getNEFoldChangePValue(Wash washPostWash, Wash washPreWash) {
		final Pair<Wash, Wash> pair = new Pair<Wash, Wash>(washPostWash, washPreWash);
		if (!nuclearEnvelopeFoldChangesPValuesPreWashVSWash.containsKey(pair)) {
			throw new IllegalArgumentException(pair + " not in the hash");
		}
		return nuclearEnvelopeFoldChangesPValuesPreWashVSWash.get(pair);
	}

	public String getProteinKey() {
		return proteinKey;
	}

	public Integer getTotalNuclearContentFractionSPCInUInSabyDataset() {
		return totalNuclearContentFractionSPCInUInSabyDataset;
	}

	public void setTotalNuclearContentFractionSPCInUInSabyDataset(
			Integer totalNuclearContentFractionSPCInUInSabyDataset) {
		this.totalNuclearContentFractionSPCInUInSabyDataset = totalNuclearContentFractionSPCInUInSabyDataset;
	}

	public Double getNuclearContentEnrichmentScoreInUInSabyDataset() {
		return nuclearContentEnrichmentScoreInUInSabyDataset;
	}

	public void setNuclearContentEnrichmentScoreInUInSabyDataset(Double nuclearContentEnrichmentScoreInUInSabyDataset) {
		this.nuclearContentEnrichmentScoreInUInSabyDataset = nuclearContentEnrichmentScoreInUInSabyDataset;
	}

}
