from os.path import isdir, join
from os import rmdir, getcwd

import pytest

from rad_tools.routines import *


class TestTerminalColours:

    def test_for_human(self):
        print("\n Please check that the following colours "
              "are displayed correctly:\n")
        print(f'{BLACK}BLACK{RESET}',
              f'{RED}RED{RESET}',
              f'{GREEN}GREEN{RESET}',
              f'{YELLOW}YELLOW{RESET}',
              f'{BLUE}BLUE{RESET}',
              f'{MAGENTA}MAGENTA{RESET}',
              f'{CYAN}CYAN{RESET}',
              f'{WHITE}WHITE{RESET}',
              sep='\n')
        assert WARNING == YELLOW
        assert OK == GREEN
        assert ERROR == RED


def test_get_256_colours():
    print("\n Please check that the following colour's table "
          "looks beautifull:\n")
    for i in range(0, 16):
        for j in range(0, 16):
            print(f'{get_256_colours(16 * i + j)}'
                  f'{16 * i + j:3}'
                  f'{RESET}',
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


def test_check_make_dir_simple():

    check_make_dir('tmp')
    assert isdir('tmp') == True

    rmdir('tmp')


def test_check_make_dir_simple_twice():

    check_make_dir('tmp')
    assert isdir('tmp') == True

    check_make_dir('tmp')
    assert isdir('tmp') == True

    rmdir('tmp')


def test_check_make_dir_complex():

    check_make_dir(join('tmp', 'tmp_inside'))
    assert isdir('tmp') == True
    assert isdir(join('tmp', 'tmp_inside')) == True

    rmdir(join('tmp', 'tmp_inside'))
    rmdir('tmp')


def test_check_make_dir_complex_twice():

    check_make_dir(join('tmp', 'tmp_inside'))
    assert isdir('tmp') == True
    assert isdir(join('tmp', 'tmp_inside')) == True

    check_make_dir(join('tmp', 'tmp_inside'))
    assert isdir('tmp') == True
    assert isdir(join('tmp', 'tmp_inside')) == True

    rmdir(join('tmp', 'tmp_inside'))
    rmdir('tmp')


def test_check_make_dir_abs_path():

    check_make_dir(join(getcwd(), 'tmp'))
    assert isdir(join(getcwd(), 'tmp')) == True

    rmdir(join(getcwd(), 'tmp'))


@pytest.mark.parametrize("iso, aniso, dmi, matrix", [
    (1,
     None,
     None,
     [[1, 0, 0],
      [0, 1., 0],
      [0, 0, 1]]),

    (1,
     [[1, 0, 0],
      [0, 1., 0],
      [0, 0, 1]],
        None,
     [[2, 0, 0],
      [0, 2., 0],
      [0, 0, 2]]),

    (1.23,
     [[0.01, 5., 0.004],
      [5., -0.003, 0],
      [0.004, 0, .034]],
        None,
        [[1.24, 5., 0.004],
         [5., 1.227, 0],
         [0.004, 0, 1.264]]),

    (None,
     [[0.01, 5., 0.004],
      [5., -0.003, 0],
      [0.004, 0, .034]],
        None,
        [[0.01, 5., 0.004],
         [5., -0.003, 0],
         [0.004, 0, 0.034]]),

    (0,
     [[0.01, 5., 0.004],
      [5., -0.003, 0],
      [0.004, 0, .034]],
        None,
        [[0.01, 5., 0.004],
         [5., -0.003, 0],
         [0.004, 0, 0.034]]),

    (None,
     [[0.01, 5., 0.004],
      [5., -0.003, 0],
      [0.004, 0, .034]],
        (0, 0, 0),
        [[0.01, 5., 0.004],
         [5., -0.003, 0],
         [0.004, 0, 0.034]]),

    (None,
     [[0.01, 5., 0.005],
      [5., -0.003, 0],
      [0.005, 0, .034]],
        (1, 0.06, -0.001),
        [[0.01, 4.999, -0.055],
         [5.001, -0.003, 1],
         [0.065, -1, 0.034]]),

    (3.0830,
     [[-0.009, 0, 0.004],
      [0, -0.013, 0],
      [0.004, 0, -0.007]],
        (-0.0002, 0.0001, 0.0001),
        [[3.074, 0.0001, 0.0039],
         [-0.0001, 3.07, -0.0002],
         [0.0041, 0.0002, 3.076]])

])
def test_exchange_to_matrix(iso, aniso, dmi, matrix):
    for i in range(0, 3):
        for j in range(0, 3):
            assert round(exchange_to_matrix(iso, aniso, dmi)[i][j], 4) ==\
                round(matrix[i][j], 4)


@pytest.mark.parametrize("iso, aniso, dmi, matrix", [
    (1,
     [[0, 0, 0],
      [0, 0, 0],
      [0, 0, 0]],
     (0, 0, 0),
     [[1, 0, 0],
      [0, 1., 0],
      [0, 0, 1]]),

    (2,
     [[0, 0, 0],
      [0, 0, 0],
      [0, 0, 0]],
        (0, 0, 0),
     [[2, 0, 0],
      [0, 2., 0],
      [0, 0, 2]]),

    (1.2437,
     [[-0.0037, 5., 0.004],
      [5., -0.0167, 0],
      [0.004, 0, 0.0203]],
        (0, 0, 0),
        [[1.24, 5., 0.004],
         [5., 1.227, 0],
         [0.004, 0, 1.264]]),

    (0.0137,
     [[-0.0037, 5., 0.004],
      [5., -0.0167, 0],
      [0.004, 0, 0.0203]],
        (0, 0, 0),
        [[0.01, 5., 0.004],
         [5., -0.003, 0],
         [0.004, 0, 0.034]]),

    (0.0137,
     [[-0.0037, 5., 0.004],
      [5., -0.0167, 0],
      [0.004, 0, 0.0203]],
        (1, 0.06, -0.001),
        [[0.01, 4.999, -0.056],
         [5.001, -0.003, 1],
         [0.064, -1, 0.034]])

])
def test_exchange_from_matrix(iso, aniso, dmi, matrix):
    iso_d, aniso_d, dmi_d = exchange_from_matrix(matrix)

    assert iso == round(iso_d, 4)

    for i in range(0, 3):
        assert round(dmi[i], 4) == round(dmi_d[i], 4)
        for j in range(0, 3):
            assert round(aniso[i][j], 4) == round(aniso_d[i][j], 4)


@pytest.mark.parametrize("mark, new_mark", [
    ('Cr1', '$Cr_{1}$'),
    ('Cr11', '$Cr_{11}$'),
])
def test_atom_mark_to_latex(mark, new_mark):
    assert atom_mark_to_latex(mark) == new_mark


@pytest.mark.parametrize("x, y, angle", [
    (1, 0, 0),
    (1, 1, 45),
    (0, 1, 90),
    (-1, 1, 135),
    (-1, 0, 180),
    (-1, -1, 225),
    (0, -1, 270),
    (1, -1, 315),
])
def test_rot_angle_dummy(x, y, angle):
    assert round(rot_angle(x, y, dummy=True), 4) == round(angle, 4)


@pytest.mark.parametrize("x, y, angle", [
    (1, 0, 0),
    (1, 1, 45),
    (0, 1, 90),
    (-1, 1, -45),
    (-1, 0, 0),
    (-1, -1, 45),
    (0, -1, 90),
    (1, -1, -45),
])
def test_rot_angle_not_dummy(x, y, angle):
    assert round(rot_angle(x, y), 4) == round(angle, 4)


def test_rot_angle_ill_case():
    with pytest.raises(ValueError):
        rot_angle(0, 0)
