from os.path import isdir, join
from os import rmdir

import pytest

from rad_tools.routines import TerminalCoulours, get_256_colours, check_make


class TestTerminalColours:

    def test_for_human(self):
        print("\n Please check that the following colours "
              "are displayed correctly:\n")
        print(f'{TerminalCoulours.BLACK}BLACK{TerminalCoulours.RESET}',
              f'{TerminalCoulours.RED}RED{TerminalCoulours.RESET}',
              f'{TerminalCoulours.GREEN}GREEN{TerminalCoulours.RESET}',
              f'{TerminalCoulours.YELLOW}YELLOW{TerminalCoulours.RESET}',
              f'{TerminalCoulours.BLUE}BLUE{TerminalCoulours.RESET}',
              f'{TerminalCoulours.MAGENTA}MAGENTA{TerminalCoulours.RESET}',
              f'{TerminalCoulours.CYAN}CYAN{TerminalCoulours.RESET}',
              f'{TerminalCoulours.WHITE}WHITE{TerminalCoulours.RESET}',
              sep='\n')
        assert TerminalCoulours.WARNING == TerminalCoulours.YELLOW
        assert TerminalCoulours.OK == TerminalCoulours.GREEN
        assert TerminalCoulours.ERROR == TerminalCoulours.RED

    def test_get_256_colours(self):
        print("\n Please check that the following colour's table "
              "looks beautifull:\n")
        for i in range(0, 16):
            for j in range(0, 16):
                print(f'{get_256_colours(16 * i + j)}'
                      f'{16 * i + j:3}'
                      f'{TerminalCoulours.RESET}',
                      end=' ')
            print()
        with pytest.raises(ValueError):
            get_256_colours(345)
        with pytest.raises(ValueError):
            get_256_colours(256)
        with pytest.raises(ValueError):
            get_256_colours(34.5)
        with pytest.raises(ValueError):
            get_256_colours(-45)
        with pytest.raises(ValueError):
            get_256_colours('sadfasd')
        assert get_256_colours(0) == '\033[38:5:0m'
        assert get_256_colours(255) == '\033[38:5:255m'
        assert get_256_colours(34) == '\033[38:5:34m'


def test_check_make():
    with pytest.raises(FileNotFoundError):
        check_make(join('sadfsd', 'sdfgsd'))
    check_make('tmp')
    assert isdir('tmp') == True
    check_make('tmp')
    rmdir('tmp')
