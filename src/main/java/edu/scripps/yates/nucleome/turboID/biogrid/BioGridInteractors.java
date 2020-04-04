package edu.scripps.yates.nucleome.turboID.biogrid;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import gnu.trove.set.hash.THashSet;

public class BioGridInteractors {

	private final Set<Interactor> interactors = new THashSet<Interactor>();
	private final List<Interactor> interactorsAsList = new ArrayList<Interactor>();
	private final String bait;

	public BioGridInteractors(String bait) {
		this.bait = bait;
	}

	public Set<Interactor> getInteractors() {
		return interactors;
	}

	public List<String> getInteractorsAsList() {
		if (interactorsAsList.isEmpty()) {
			interactorsAsList.addAll(interactors);
		}
		return interactorsAsList.stream().map(interactor -> interactor.getId().toUpperCase())
				.collect(Collectors.toList());
	}

	public boolean addInteractor(String interactor, Double score1, Double score2) {
		return this.interactors.add(new Interactor(interactor, score1, score2));
	}

	public String getBait() {
		return bait;
	}

	public Set<String> getInteractorByMinScore(double score1Threshold, double score2Threshold) {
		getInteractorsAsList();

		Collections.sort(interactorsAsList, new Comparator<Interactor>() {

			@Override
			public int compare(Interactor o1, Interactor o2) {
				return o2.getScore1().compareTo(o1.getScore1());
			}
		});
		// the first one will have the biggest score
		int max = interactorsAsList.size() - 1;
		int min = 0;
		int index = findIndexInFirstScore(interactorsAsList, score1Threshold, min, max);
		List<Interactor> subList = interactorsAsList.subList(0, index + 1);
		// now the second score
		// first remove the ones with score2=NaN
		subList = subList.stream().filter(i -> !Double.isNaN(i.getScore2())).collect(Collectors.toList());
		Collections.sort(subList, new Comparator<Interactor>() {

			@Override
			public int compare(Interactor o1, Interactor o2) {
				if (Double.isNaN(o1.getScore2())) {
					return 1;
				}
				if (Double.isNaN(o2.getScore2())) {
					return -1;
				}
				return o1.getScore2().compareTo(o2.getScore2());
			}
		});
		// the first one will have the lowest score
		max = subList.size() - 1;
		min = 0;
		if (!subList.isEmpty()) {
			index = findIndexInSecondScore(subList, score2Threshold, min, max);
			subList = subList.subList(0, index + 1);
		}
		final Set<String> set = new THashSet<String>();

		for (final Interactor interactor : subList) {
			set.add(interactor.getId().toUpperCase());
		}
		return set;
	}

	private int findIndexInFirstScore(List<Interactor> list, double scoreThreshold, int min, int max) {
		if (min == max) {
			return min;
		}
		if (min == max - 1) {
			return min;
		}
		final int index = (max - min) / 2 + min;
//		System.out.println(min + "\t" + max + "\t" + index + "\t" + list.get(index).getScore());

		if (list.get(index).getScore1() < scoreThreshold) {
			// look into the first half
			return findIndexInFirstScore(list, scoreThreshold, min, index);
		} else {
			// look into the second half
			return findIndexInFirstScore(list, scoreThreshold, index, max);
		}
	}

	private int findIndexInSecondScore(List<Interactor> list, double scoreThreshold, int min, int max) {
		if (min == max) {
			return min;
		}
		if (min == max - 1) {
			return min;
		}
		final int index = (max - min) / 2 + min;
//		System.out.println(min + "\t" + max + "\t" + index + "\t" + list.get(index).getScore());

		if (list.get(index).getScore2() > scoreThreshold) {
			// look into the first half
			return findIndexInSecondScore(list, scoreThreshold, min, index);
		} else {
			// look into the second half
			return findIndexInSecondScore(list, scoreThreshold, index, max);
		}
	}
}
