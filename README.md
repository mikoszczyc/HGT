# HGT pseudocode 

    for seq_id in alignment.identity:

	# seq_id = species1|protein1
	# other_seq -> sequences other than seq_id (from alignment.identity file)
	identity = [(other_seq1, %identity1), (other_seq2, %identity2), ...]
	sort identity

	seq_id.genus = ...
	seq_id.family = ...
	seq_id.order = ...

	found = false
	for el in identity:
		if el ∉ genus ∧ el ∈ family:
			found = true
			hit = el / hit = (el_seq, %el_identity)
			break

	if not found:
		for el in identity:
			if el ∉ family ∧ el ∈ order:
				found = true
				hit = el / hit = (el_seq, %el_identity)
				break

    if not found:
		for el in identity:
			if el ∉ order ∧ el ∈ class:
				found = true
				hit = el / hit = (el_seq, %el_identity)
				break

	save all hits to 'hit.txt'

	sort hits, pair them and save pairs to 'crossedResult.csv'