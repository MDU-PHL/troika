import sys, pathlib, pandas, pytest, numpy

from unittest.mock import patch

from troika.RunTroika import RunTroika



def test_3_cols():
        '''
        return True when correct number of columns
        '''
        with patch.object(RunTroika, "__init__", lambda x: None):
                detect_obj = RunTroika()
                tab = pandas.DataFrame({'A':[1], 'B':[2], 'C':[3]})
                assert detect_obj.three_cols(tab)


def test_2col_dimensions():
        '''
        return False when wrong number of columns present
        '''
        with patch.object(RunTroika, "__init__", lambda x: None):
                detect_obj = RunTroika()
                tab = pandas.DataFrame({'A':[1], 'B':[2]})
                assert detect_obj.three_cols(tab) == False


def test_path_exists():
        '''
        test that path_exists returns True
        '''
        with patch.object(RunTroika, "__init__", lambda x: None):
                p = pathlib.Path('troika','tests', 'troika.txt')
                detect_obj = RunTroika()
                assert detect_obj.path_exists(p)

