package edu.scripps.yates.nucleome.turboID.biogrid;

import org.apache.commons.lang3.builder.HashCodeBuilder;

public class Interactor {
	private final String id;
	private final Double score1;
	private final Double score2;
	private Integer hashCode;

	public Interactor(String id, Double score1, Double score2) {
		this.id = id;
		this.score1 = score1;
		this.score2 = score2;
	}

	public String getId() {
		return id;
	}

	public Double getScore1() {
		return score1;
	}

	public Double getScore2() {
		return score2;
	}

	@Override
	public String toString() {
		return id;
	}

	@Override
	public int hashCode() {
		if (hashCode == null) {
			hashCode = HashCodeBuilder.reflectionHashCode(id, false);
		}
		return hashCode;
	}

	@Override
	public boolean equals(Object obj) {
		if (obj instanceof Interactor) {
			final Interactor interactor = (Interactor) obj;
			if (interactor.getId().equalsIgnoreCase(getId())) {
				return true;
			}
			return false;
		}
		return super.equals(obj);
	}
}
