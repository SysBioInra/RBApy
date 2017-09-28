"""Module defining CurationData class."""

# python 2/3 compatibility
from __future__ import division, print_function

import pandas


class CurationData(object):
    """
    Class used to read/write data in a csv file.

    Attributes:
        missing_tag: tag used when some piece of information is missing.
        data: pandas.DataFrame holding current information.

    """

    missing_tag = '[MISSING]'

    def __init__(self, columns):
        """
        Constructor.

        Args:
            columns: list of default column names.
        """
        self.data = pandas.DataFrame(columns=columns)

    def read(self, input_file):
        """
        Read data from file, overriding existing data.

        Args:
            input_file: path to file or buffer.
        """
        self.data = pandas.read_csv(input_file, sep='\t',
                                    na_values=[self.missing_tag],
                                    keep_default_na=False)

    def write(self, output_file):
        """
        Write data to file.

        Args:
            output_file: path to file or buffer.
        """
        self.data.to_csv(output_file, sep='\t', na_rep=self.missing_tag,
                         index=False)

    def add(self, rows):
        """
        Add rows to data.

        Args:
            rows: set of rows to add to current data. The length of each row
                has to match the number of columns.
        """
        new_rows = pandas.DataFrame(rows, columns=self.data.columns)
        self.data = self.data.append(new_rows)

    def has_missing_information(self):
        """
        Check whether any piece of information is missing.

        Returns:
            True if any field has missing information.

        """
        return pandas.isnull(self.data).values.any()
