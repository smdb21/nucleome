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

	public boolean addInteractor(String interactor, Double score) {
		return this.interactors.add(new Interactor(interactor, score));
	}

	public String getBait() {
		return bait;
	}

	public Set<String> getInteractorByMinScore(double scoreThreshold) {
		getInteractorsAsList();

		Collections.sort(interactorsAsList, new Comparator<Interactor>() {

			@Override
			public int compare(Interactor o1, Interactor o2) {
				return o2.getScore().compareTo(o1.getScore());
			}
		});
		// the first one will have the biggest score
		final int max = interactorsAsList.size() - 1;
		final int min = 0;
		final int index = findIndex(interactorsAsList, scoreThreshold, min, max);
		final Set<String> set = new THashSet<String>();
		final List<Interactor> subList = interactorsAsList.subList(0, index + 1);
		for (final Interactor interactor : subList) {
			set.add(interactor.getId().toUpperCase());
		}
		return set;
	}

	private int findIndex(List<Interactor> list, double scoreThreshold, int min, int max) {
		if (min == max) {
			return min;
		}
		if (min == max - 1) {
			return min;
		}
		final int index = (max - min) / 2 + min;
//		System.out.println(min + "\t" + max + "\t" + index + "\t" + list.get(index).getScore());

		if (list.get(index).getScore() < scoreThreshold) {
			// look into the first half
			return findIndex(list, scoreThreshold, min, index);
		} else {
			// look into the second half
			return findIndex(list, scoreThreshold, index, max);
		}
	}
}
