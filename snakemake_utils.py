import re


def config_key_or_fail(config, key):
	"""Get the value of a config key, raise informative error if does not exist.

	This is a bit better than just raising a KeyError when using basic __getitem__ syntax.
	"""
	try:
		return config[key]
	except KeyError:
		raise RuntimeError('Missing config key {!r}'.format(key)) from None


def find_read_file_pairs(files):
	"""Find pairs of read files in a set of file names.

	Based on Tom de Man's scripts.
	"""

	pairs = dict()
	unpaired = []

	# Find pairs
	for file in files:

		match = re.match(r'(.+)_R?(1|2)[_.$].*', file)

		if match is not None:

			name = match.group(1)
			index = int(match.group(2)) - 1

			pair = pairs.setdefault(name, [None, None])

			if pair[index] is not None:
				raise ValueError(
					'Member of file pair identified twice: {}, {}'
					.format(file, pair[index])
				)

			pair[index] = file

		else:
			unpaired.append(file)

	# Throw out pairs with only one item found
	for key, pair in list(pairs.items()):
		if None in pair:
			del pairs[key]
			for file in pair:
				if file is not None:
					unpaired.append(file)

	return {k: tuple(v) for k, v in pairs.items()}, unpaired



def parse_fasta_basic(lines):

	current_id = None
	current_seq = None

	for line in lines:
		line = line.strip()

		# Skip blank lines and comments
		if not line or line.startswith('#'):
			continue

		if line.startswith('>'):
			# Start of new sequence

			if current_id is not None:
				yield current_id, current_seq

			current_id = line[1:]
			current_seq = ''

		else:
			# Assume sequence data
			if current_id is None:
				raise ValueError('Sequence data before first ID line')

			current_seq += line

	if current_id is not None:
		yield current_id, current_seq


def write_fasta_seq(fobj, id_, seq, wrap=80):

	# Write ID line
	fobj.write('>' + id_ + '\n')

	# Write lines one at a time up to line width
	for i in range(0, len(seq), wrap):
		fobj.write(seq[i:i + wrap] + '\n')
