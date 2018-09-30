"""Module defining CurationData class."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

import pandas


class CurationData(object):
    """
    Read/write data in a csv file.

    Parameters
    ----------
    missing_tag : str
        tag used when some piece of information is missing.
    data : pandas.DataFrame
        data.
    filename : str
        Path to file used to read/write data

    """

    missing_tag = '[MISSING]'

    def __init__(self, filename, columns):
        """
        Build object from existing curation file or create one.

        Parameters
        ----------
        filename : str
            Location of curation file.
        columns : list of str
            Names of columns of curation file. Names will be overridden if
            a curation file already exists.

        """
        self.filename = filename
        self._data_added = False
        self.data = pandas.DataFrame(columns=columns)
        try:
            self.data = pandas.read_csv(filename, sep='\t',
                                        na_values=[self.missing_tag])
        except IOError:
            print('Helper file {} not found.'.format(filename))
            self.write(self.filename)

    def update_file(self):
        """
        Write data to default file if data was added since last update.

        Returns
        -------
        bool
            True if file needed to be updated.

        """
        if self._data_added:
            self.write(self.filename)
            self._data_added = False
            return True
        return False

    def write(self, output_file):
        """
        Write data to file.

        Parameters
        ----------
        output_file : str
            Path to file or buffer.

        """
        self.data.to_csv(output_file, sep='\t',
                         na_rep=self.missing_tag, index=False)

    def rows(self):
        """
        Access to data as a simple array.

        Returns
        -------
        numpy.array
            Array containing data.

        """
        return self.data.values

    def add_rows(self, rows):
        """
        Add rows to data.

        Parameters
        ----------
        rows : list of lists/tuples
            Set of rows to add to current data. The length of each row
            has to match the number of columns.

        """
        new_rows = pandas.DataFrame(rows, columns=self.data.columns)
        self.data = self.data.append(new_rows)
        self._data_added = True

    def add_row(self, row):
        """
        Add single row to data.

        Parameters
        ----------
        row : list/tuples
            Row to add to current data.

        """
        self.add_rows([row])

    def has_missing_information(self, columns=None):
        """
        Check whether any piece of information is missing.

        Parameters
        ----------
        columns : str or list of str
            Columns to check. If None (the default), check all columns.

        Returns
        -------
        bool
            True if any column has missing information.

        """
        if columns:
            return pandas.isnull(self.data[columns]).values.any()
        else:
            return pandas.isnull(self.data).values.any()
