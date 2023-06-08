"""A module containing classes for aligning structure."""

import copy
import math
import random
import sys
from typing import List, Optional

import numpy

from .geometry import unit_vector
from .protein import Polypeptide


def align_backbones(reference, mobile, stop_when=None, verbose=False):
    mobile = copy.deepcopy(mobile)
    initial_trans = reference.centre_of_mass - mobile.centre_of_mass
    mobile.translate(initial_trans)
    fitter = MMCAlign(_align_eval, [reference], mobile)
    fitter.start_optimisation(
        500, 10, 1, temp=100, stop_when=stop_when, verbose=verbose
    )
    return fitter.best_energy


def _align_eval(loop, reference):
    return loop.rmsd(reference, backbone=True)


class MMCAlign:
    """A alignment protocol that uses Metropolis Monte Carlo.

    Notes
    -----
    THIS IS CURRENTLY SUPER INEFFICIENT DUE TO THE DEEPCOPIES.
    I plan to improve this by aligning arrays of atoms and only
    recording the rotation and translation that led to that alignment.

    Parameters
    ----------
    eval_fn : Polypeptide -> float
        A function to evaluate the quality of your fit.
    eval_args : list
        A list of static args to be used in the `eval_fn`, these will
        be unpacked into the evaluation function _i.e._
        `eval_fn(polypeptide, *eval_args).
    polypeptide : Polypeptide
        An ampal polypeptide containing the model to be aligned.
    """

    def __init__(
        self, eval_fn, eval_args: Optional[list], polypeptide: Polypeptide
    ) -> None:
        self.eval_fn = eval_fn
        if eval_args is None:
            self.eval_args: List = []
        else:
            self.eval_args = eval_args
        self.current_energy = None
        self.best_energy = None
        self.best_model = None
        self.polypeptide = polypeptide

    def start_optimisation(
        self,
        rounds: int,
        max_angle: float,
        max_distance: float,
        temp: float = 298.15,
        stop_when=None,
        verbose=None,
    ):
        """Starts the loop fitting protocol.

        Parameters
        ----------
        rounds : int
            The number of Monte Carlo moves to be evaluated.
        max_angle : float
            The maximum variation in rotation that can moved per
            step.
        max_distance : float
            The maximum distance the can be moved per step.
        temp : float, optional
            Temperature used during fitting process.
        stop_when : float, optional
            Stops fitting when energy is less than or equal to this value.
        """
        self._generate_initial_score()
        self._mmc_loop(
            rounds,
            max_angle,
            max_distance,
            temp=temp,
            stop_when=stop_when,
            verbose=verbose,
        )
        return

    def _generate_initial_score(self):
        """Runs the evaluation function for the initial pose."""
        self.current_energy = self.eval_fn(self.polypeptide, *self.eval_args)
        self.best_energy = copy.deepcopy(self.current_energy)
        self.best_model = copy.deepcopy(self.polypeptide)
        return

    def _mmc_loop(
        self, rounds, max_angle, max_distance, temp=298.15, stop_when=None, verbose=True
    ):
        """The main Metropolis Monte Carlo loop."""
        current_round = 0
        while current_round < rounds:
            working_model = copy.deepcopy(self.polypeptide)
            random_vector = unit_vector(numpy.random.uniform(-1, 1, size=3))
            mode = random.choice(["rotate", "rotate", "rotate", "translate"])
            if mode == "rotate":
                random_angle = numpy.random.rand() * max_angle
                working_model.rotate(
                    random_angle, random_vector, working_model.centre_of_mass
                )
            else:
                random_translation = random_vector * (
                    numpy.random.rand() * max_distance
                )
                working_model.translate(random_translation)
            proposed_energy = self.eval_fn(working_model, *self.eval_args)
            move_accepted = self.check_move(
                proposed_energy, self.current_energy, t=temp
            )
            if move_accepted:
                self.current_energy = proposed_energy
                if self.current_energy < self.best_energy:
                    self.polypeptide = working_model
                    self.best_energy = copy.deepcopy(self.current_energy)
                    self.best_model = copy.deepcopy(working_model)
            if verbose:
                sys.stdout.write(
                    "\rRound: {}, Current RMSD: {}, Proposed RMSD: {} "
                    "(best {}), {}.       ".format(
                        current_round,
                        self.float_f(self.current_energy),
                        self.float_f(proposed_energy),
                        self.float_f(self.best_energy),
                        "ACCEPTED" if move_accepted else "DECLINED",
                    )
                )
                sys.stdout.flush()
            current_round += 1
            if stop_when:
                if self.best_energy <= stop_when:
                    break
        return

    @staticmethod
    def float_f(f):
        """Formats a float for printing to std out."""
        return "{:.2f}".format(f)

    @staticmethod
    def check_move(new, old, t):
        """Determines if a model will be accepted."""
        if (t <= 0) or numpy.isclose(t, 0.0):
            return False
        K_BOLTZ = 1.9872041e-003  # kcal/mol.K
        if new < old:
            return True
        else:
            move_prob = math.exp(-(new - old) / (K_BOLTZ * t))
            if move_prob > random.uniform(0, 1):
                return True
        return False


__author__ = "Christopher W. Wood"
