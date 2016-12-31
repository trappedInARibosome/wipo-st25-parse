# Copyright 2006-2015 by Chris Jackson.  All rights reserved.
# Please see the LICENSE file that should have been included
# as part of this package.
#
# This module is for reading WIPO ST.25 format files as SeqRecord objects.
# WIPO ST.25 (USPTO 2422) standard files are designed to be human readable,
# and are frequently not fully compliant with the published specification.
#
# This package is provided entirely as-is, and should work if your file meets
# the published specification. There is no guarantee or warranty.

from __future__ import print_function

from Bio.Alphabet import generic_protein
from Bio.Alphabet import generic_nucleotide
from Bio.Alphabet import generic_rna

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation

from Bio.SeqIO.Interfaces import SequenceIterator

# Protein three to one letter conversion table
THREE_LETTER_PROT_DICT = {"ala": "A",
                          "arg": "R",
                          "asn": "N",
                          "asp": "D",
                          "cys": "C",
                          "glu": "E",
                          "gln": "Q",
                          "gly": "G",
                          "his": "H",
                          "ile": "I",
                          "leu": "L",
                          "lys": "K",
                          "met": "M",
                          "phe": "F",
                          "pro": "P",
                          "ser": "S",
                          "thr": "T",
                          "trp": "W",
                          "tyr": "Y",
                          "val": "V",
                          "xaa": "X",
                          "glx": "Z",
                          "asx": "B"}

ONE_LETTER_NUCL_LIST = list("agcturymkswbdhvn")

# Maps these numeric identifier field codes to the corresponding seqrecord object attribute
SEQRECORD_MAP = {"<210>": "id",
                 "<213>": "description"}

# Map these numeric identifier field codes to seqrecord object attributes
HEADER_MAP = {"<110>": "applicant",
              "<120>": "title",
              "<130>": "file_reference",
              "<140>": "application_number",
              "<141>": "file_date",
              "<150>": "prior_application",
              "<160>": "seq_count",
              "<170>": "software"}

# Attribute field to consider the start of a new record
RECORD_START_KEY = "<210>"

# Attribute field to parse for sequence type
SEQUENCE_TYPE_KEY = "<212>"

# Attribute field to parse for sequence
SEQUENCE_KEY = "<400>"

# Attribute field to parse for sequence length
SEQUENCE_LENGTH_KEY = "<211>"

# Map these sequence types to the corresponding Bio.Alphabets
SEQUENCE_TYPE_MAP = {"PRT": generic_protein,
                     "DNA": generic_nucleotide,
                     "RNA": generic_rna}

# Attribute field to consider the start of a new feature set
FEATURE_KEY = "<220>"

# Map these numeric identifier field codes to seqfeature object attributes
FEATURE_MAP = {"<221>": "type",
               "<223>": "other_information"}

FEATURE_LOCATION_KEY = "<222>"


class ST25SequenceIterator(SequenceIterator):
    def __init__(self, handle, check_length=True, seq_len_key=SEQUENCE_LENGTH_KEY, seq_key=SEQUENCE_KEY,
                 header_map=HEADER_MAP, seqrecord_map=SEQRECORD_MAP, feature_key=FEATURE_KEY,
                 translation_table=THREE_LETTER_PROT_DICT, allowed_nucl_char=ONE_LETTER_NUCL_LIST,
                 feature_map=FEATURE_MAP, location_key=FEATURE_LOCATION_KEY, type_key=SEQUENCE_TYPE_KEY,
                 type_map=SEQUENCE_TYPE_MAP, break_key=RECORD_START_KEY):

        """
        ST25SequenceIterator takes the ST25 file handle, and any of the file constants that you'd want to overwrite, and
        acts as an generator for SeqRecord objects.

        This whole thing is kinda a really hacked together mess, but so is the file format so we're even

        Required (Positional) Arguments:

        :param handle: handle
            An open file handle to the target file

        Optional (Keyword) Arguments:

        :param seq_len_key: str
            Numeric identifier for the sequence length field. Defaults to SEQUENCE_LENGTH_KEY

        :param seq_key: str
            Numeric identifier for the sequence field. Defaults to SEQUENCE_KEY

        :param header_map: dict
            Dict to map numeric identifiers in the header block to object attributes. Defaults to HEADER_MAP

        :param seqrecord_map: dict
            Dict to map numeric identifiers in the sequence block to object attributes. Defaults to SEQRECORD_MAP

        :param feature_key: str
            Numeric identifier for the start of a new feature block. Defaults to FEATURE_KEY

        :param header_map: dict
            Map of the numeric identifiers to SeqRecord object attributes. Defaults to HEADER_MAP

        :param translation_table: dict
            A dict for translating three letter protein codes to one letter protein sequence. Keys should be all lower
            case. Defaults to THREE_LETTER_PROT_DICT

        :param allowed_nucl_char: list
            A list of characters to keep as nucleotide sequence. Defaults to ONE_LETTER_NUCL_LIST

        :param feature_map: dict
            A dict for mapping feature numeric identifiers to SeqFeature object attributes. Defaults to FEATURE_MAP

        :param location_key: str
            The numeric identifier for location fields. Defaults to FEATURE_LOCATION_KEY

        :param type_key: str
            The key that will be used to identify the sequence type record, so that the sequence alphabet can be
            determined. Defaults to the SEQUENCE_TYPE_KEY constant

        :param type_map: dict
            A dict for mapping sequence types to Bio.Alphabet objects

        :param break_key: str
            The numeric identifier to use as the indication that a new record has started. Defaults to BREAK_RECORD_KEY
        """

        SequenceIterator.__init__(self, handle)

        # Time to write some constants/keyword arguments
        self.check_length = check_length
        self.seq_len_key = seq_len_key
        self.seq_key = seq_key
        self.header_map = header_map
        self.seqrecord_map = seqrecord_map
        self.feature_key = feature_key
        self.translation_table = translation_table
        self.allowed_nucl_char = allowed_nucl_char
        self.feature_map = feature_map
        self.location_key = location_key
        self.type_key = type_key
        self.type_map = type_map
        self.break_key = break_key

        self.header_data = None

    def __iter__(self):

        """
        Return self as an iterator object

        :return: ST25SequenceIterator
        """

        self.chunk_generator = self._chunk_records()
        return self

    def __next__(self):

        """
        Keeps looping until it gets a usable SeqRecord object or until a StopIteration is raised from chunk_generator

        :return: SeqRecord
        """

        while True:

            # Get a file chunk and then try parsing it
            file_chunk = next(self.chunk_generator)

            try:
                # If parsing is successful, break and return. Other errors can be raised which are not caught here.
                seq_record = self.parse_st25(file_chunk)
                break

            # Roll again if there's a TypeError (this is the error raised for non-sequence or non-parsable file chunks)
            except TypeError:
                continue

            # Reraise StopIteration when the chunk_generator hits the end of the file
            except StopIteration:
                raise

        return seq_record

    def parse_st25(self, file_chunk):
        """
        parse_st25 is the main worker function to take a file chunk, which is a list of lines as strings, and turn it
        into a Bio.SeqRecord object. It will raise a TypeError if the file chunk doesnt look like sequence, and will
        raise a ValueError if there is a mismatch between the sequence length in the file and the parsed sequence length

        This whole thing is kinda a really hacked together mess, but so is the file format so we're even

        Required (Positional) Arguments:

        :param file_chunk: list [str]
            A chunk of the input file handle as a list of strings

        Returns:

        :return: Bio.SeqRecord
            Returns a SeqRecord object
        """

        record_chunk, chunk_type = self._st25_process_chunk(file_chunk)

        # If this file chunk has the same keys as HEADER_DICT and the header_data hasn't been set, parse it as a
        # header and then set header_data. Yea this is super hacky what can you do.
        if self.header_data is None and chunk_type is None:
            sum_header_lines = sum(list(record_chunk.keys()).count(head_i) for head_i in self.header_map.keys())
            if sum_header_lines > 1:
                self.header_data = self._st25_parse_header(record_chunk)
                raise TypeError("Header block")

        # Use the alphabet object in chunk_type to pass the sequence block to protein or nucleotide parser
        if chunk_type is generic_protein:
            seq_obj = self._st25_parse_protein(record_chunk[self.seq_key], alphabet=chunk_type)
        elif chunk_type is generic_nucleotide or chunk_type is generic_rna:
            seq_obj = self._st25_parse_nucl(record_chunk[self.seq_key], alphabet=chunk_type)
        else:
            raise TypeError("Unknown sequence")

        seq_record = SeqRecord(seq_obj)

        # If there's header data, use the header map to set SeqRecord object attributes with the appropriate data
        if self.header_data is not None:
            for name_id in self.header_data:
                setattr(seq_record, name_id, self.header_data[name_id])

        # If there's data fields with numeric identifiers in the seqrecord_map, set SeqRecord object attributes with the
        # appropriate data
        for numeric_id in self.seqrecord_map.keys():
            try:
                attr_to_set = ""
                for record_line in record_chunk[numeric_id]:
                    attr_to_set += " " + self._st25_get_value(record_line)
                setattr(seq_record, self.seqrecord_map[numeric_id], attr_to_set.strip())
            except KeyError:
                pass

        # If there are any features, turn them into SeqFeature objects with FeatureLocation locations (if present)
        features = []
        try:
            for feature_dict in record_chunk[self.feature_key]:
                features.append(self._st25_parse_feature(feature_dict))
        except KeyError:
            pass
        if len(features) > 0:
            setattr(seq_record, "features", features)

        # If the check_length flag is set, make sure that the length in the file matches the parsed sequence object len
        if self.check_length:
            try:
                seq_quoted_len = int(self._st25_get_value(record_chunk[self.seq_len_key][0]))

                if seq_quoted_len != len(seq_obj):
                    raise ValueError(
                        "Sequence length mismatch in ID: {} (File annotation length {}, Parsed length {}".format(
                            seq_record.id,
                            seq_quoted_len,
                            len(seq_obj)))
            except IndexError:
                pass

        return seq_record

    def _st25_parse_header(self, header_lines):
        """
        Parse header block

        Required (Positional) Arguments:

        :param header_lines: dict
            File lines keyed by numeric identifier from _st25_process_chunk

        Returns:

        :return: dict
            Header data values keyed by the seqrecord attribute name that they should populate
        """

        header_data = {}

        for numeric_id, file_lines in header_lines.items():

            values = []

            for line in file_lines:
                values.append(self._st25_get_value(line))

            try:
                header_data[self.header_map[numeric_id]] = " ".join(values)
            except KeyError:
                raise ValueError("Unknown numeric identifier '{}' in header".format(numeric_id))

        return header_data

    def _st25_parse_protein(self, sequence_lines, alphabet=generic_protein):
        """
        Parse protein-type data into a Seq object. Splits each line by whitespace and then looks for things that are in
        the translation table key list. Throws away anything else. Also throws away any line with a colon in it just
        to be on the safe side

        Required (Positional) Argument:

        :param sequence_lines: list
            List of file lines as strings from a single protein sequence record

        Optional (Keyword) Arguments:

        :param alphabet: Bio.Alphabet

        Returns:

        :return: Bio.Seq
            Seq object with the parsed single-letter protein sequence
        """

        active_sequence_string = ""

        for line in sequence_lines:

            if ":" in line:
                continue

            line_tabs = line.strip().split()

            if len(line) < 2 or len(line_tabs) < 1:
                continue

            for aa_triple_code in line_tabs:
                try:
                    active_sequence_string += self.translation_table[aa_triple_code.lower()]
                except KeyError:
                    pass

        return Seq(active_sequence_string, alphabet=alphabet)

    def _st25_parse_nucl(self, sequence_lines, alphabet=generic_nucleotide):
        """
        Parse nucleotide data into a Seq object. Throws away any line with a colon in it just to be on the safe side,
        then keeps any character in the allowed character list

        Required (Positional) Argument:

        :param sequence_lines: list
            List of file lines as strings from a single nucleotide sequence record

        Optional (Keyword) Arguments:

        :param alphabet: Bio.Alphabet

        Returns:

        :return: Bio.Seq
            Seq object with the parsed nucleotide sequence
        """

        active_sequence_string = ""

        for line in sequence_lines:

            if ":" in line:
                continue

            for char in line.strip():
                if char.lower() in self.allowed_nucl_char:
                    active_sequence_string += char

        return Seq(active_sequence_string, alphabet=alphabet)

    def _st25_parse_feature(self, feature_dict):
        """
        Parse a block of feature data into a SeqFeature object. If there's a location block, turn it into a
        FeatureLocation object.

        Required (Positional) Argument:

        :param feature_dict: dict
            Dict of numeric identifiers that are keying a list of file lines

        Returns:

        :return: Bio.SeqFeature
        """

        feature = SeqFeature()

        # Take all the keys in feature_map and, for any of them that exist in feature_dict, get the field values and
        # set them as attributes in the SeqFeature object
        for feat_key in self.feature_map.keys():
            try:
                attr_to_set = " ".join(map(self._st25_get_value, feature_dict[feat_key]))
                setattr(feature, self.feature_map[feat_key], attr_to_set)
            except KeyError:
                pass

        # If the location field is set, pass the location information to _st25_parse_location and then set the location
        # attribute in the SeqFeature object with a list of FeatureLocations
        locations = []
        try:
            for loc in feature_dict[self.location_key]:
                seq_loc = self._st25_parse_location(loc)
                if seq_loc is not None:
                    locations.append(seq_loc)
        except KeyError:
            pass

        if len(locations) > 0:
            setattr(feature, "location", locations)

        return feature

    def _st25_parse_location(self, location_line, location_middle=".."):
        """
        Parse a location line into a FeatureLocation .

        Required (Positional) Argument:

        :param location_line: str
            The file line with the location numeric identifier

        Optional (Keyword) ArgumentsL

        :param location_middle: str
            The string to break the location line at

        Returns:

        :return: Bio.SeqFeature.FeatureLocation
        """

        assert (isinstance(location_line, str))

        locations = self._st25_get_value(location_line).split(location_middle)

        if len(locations) != 2:
            return None

        start = int(_strip_non_digit(locations[0]))
        stop = int(_strip_non_digit(locations[1]))

        # Get the strand and make sure that the lowest location number is set as the feature start
        if start > stop:
            seq_loc = FeatureLocation(stop, start, strand=-1)
        else:
            if start == stop:
                strand = 0
            else:
                strand = 1
            seq_loc = FeatureLocation(start, stop, strand=strand)

        return seq_loc

    def _st25_process_chunk(self, chunk_list_of_strings):
        """
        _st25_process_chunk takes a list of file lines that comprise a single record and returns it as a dict keyed by
        line numeric identifier, as defined by WIPO ST.25 (See USPTO 2422). It processes the SEQUENCE_TYPE_KEY record
        into a Bio.Alphabet, and it processes any feature records into a list of feature dicts, each of which is keyed
        by the feature numeric identifier

        Positional (Required) Arguments:

        :param chunk_list_of_strings: list [str]
            A list of file lines for a sequence block

        Returns:

        :return processed_data_dict: dict
            A dict keyed by numeric identifier (including angle brackets, so <210>) containing a list with individual
            lines from that numeric identifier as strings.

        :return sequence_alphabet: Bio.Alphabet

        """

        feature_keys = list(self.feature_map.keys()) + [self.location_key]

        # Holds multiline numeric identifiers during scan
        active_key = None
        active_feature = None

        # Holds the data to be returned
        processed_data_dict = {}
        sequence_alphabet = None

        for line in chunk_list_of_strings:

            line_key = self._st25_get_numeric_id(line)

            # If there's no line numeric identifier, and none has been seen yet, just keep moving
            if line_key is None and active_key is None:
                continue

            # If there's a numeric identifier, set it to active key and save the line contents for parsing
            if line_key is not None:
                active_key = line_key
                line_contents = line.strip()[5:]
            else:
                line_contents = line.strip()

            # If the numeric identifier matches the sequence type identifier, parse the line and set the sequence
            # alphabet
            if active_key == self.type_key:
                try:
                    sequence_alphabet = self.type_map[self._st25_get_value(line_contents)]
                except KeyError:
                    sequence_alphabet = None

            # If the numeric identifier matches the feature identifier,
            elif active_key == self.feature_key:
                if active_feature is not None:
                    _add_dict_list(processed_data_dict, active_key, active_feature)

                active_feature = {active_key: line_contents}

            elif active_key in feature_keys:
                _add_dict_list(active_feature, active_key, line_contents)

            else:
                _add_dict_list(processed_data_dict, active_key, line_contents)

        if active_feature is not None:
            _add_dict_list(processed_data_dict, self.feature_key, active_feature)

        return processed_data_dict, sequence_alphabet

    def _chunk_records(self):
        """
        _chunk_records takes the file handle and yields a list of lines broken at a specified numeric identifier key.
        By default it breaks the records on the <210> SEQ ID NO line.

        Yields:

        :yield active_chunk: list [str]
            A list of file lines from the handle, broken at break_key
        """

        active_chunk = []

        for line in self.handle:

            line_key = self._st25_get_numeric_id(line)

            # If the line numeric identifier is the ID to break the records at, then yield the most recent chunk
            # And then start a new chunk
            if line_key == self.break_key:

                return_chunk = active_chunk
                active_chunk = [line.strip()]

                yield return_chunk

            else:

                active_chunk.append(line.strip())

        yield active_chunk

    def _st25_get_numeric_id(self, line_to_check):
        assert isinstance(line_to_check, str)

        try:
            line_key = line_to_check.strip().split()[0]
        except IndexError:
            return None

        # Hacky workaround to some files not having whitespace
        if line_key[0] == "<" and len(line_key) > 5:
            line_key = line_to_check.strip()[:5]

        if line_key[0] != "<":
            return None
        else:
            return line_key

    def _st25_get_value(self, line_to_check):
        assert isinstance(line_to_check, str)

        line_split = line_to_check.strip().split(":")
        if len(line_split) > 1:
            return line_split[-1].strip()
        else:
            return line_to_check.strip().split()[-1]


def _add_dict_list(dict_to_add, key, value):
    try:
        dict_to_add[key].append(value)
    except KeyError:
        dict_to_add[key] = [value]


def _strip_non_digit(string):
    return "".join(char for char in string if char.isdigit())
