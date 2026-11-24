"""Shape work at fixed active labels with the R4 melted labels in harmonic equilibrium."""
import numpy as np

from mouth_cp10_field_shape import Bath
from mouth_free_shape_solid import Solid
from mouth_r4_eulerian_relax_labels import InfeasibleLabels, LabelRelaxation, chart_orientation
from mouth_r4_eulerian_release import Driver, block_axes, order_grid
from mouth_r4_eulerian_solid import InvalidSpatialCrystal, SpatialCrystal
from mouth_r4_fold_order import axes_model


class ReducedShapeWork:
    def __init__(self, driver, trial, solid):
        self.driver = driver
        self.trial = trial
        self.solid = solid
        self.initial_shape = driver.x.copy()
        self.order = order_grid(trial,driver.fields)
        self.bath = Bath(trial,driver.base,driver.fields,driver.reference_depth)
        active = LabelRelaxation(trial.b,driver.labels,self.order,driver.epsilon)
        self.trace = active._pack(driver.labels)
        self.last_dx = self.last_labels = self.last_info = None

    def __call__(self, dx):
        driver, trial = self.driver, self.trial
        if np.any(dx):
            shape = driver.shape(self.initial_shape + dx)
            body = Solid(driver.prep,shape,trial.polar_shift+dx[5],8,24,6,
                         axes_model(shape,driver.reference_depth),degree=3,subdivisions=2)
            solid = SpatialCrystal(body)
        else:
            body, solid = trial.b, self.solid
        # Free trace entries are offsets from the normalized identity chart.
        # Holding the packed vector fixed implements transport onto the new knots.
        harmonic = LabelRelaxation(body,driver.labels,self.order,driver.epsilon)
        harmonic._components()
        labels = harmonic._lift(self.trace)
        crystal, _, _, orientation = solid.evaluate(labels,self.order,False)
        if orientation['min_chart_Jinv_bound'] < np.sqrt(np.finfo(float).eps):
            try:
                chart_orientation(block_axes(body),labels,harmonic.active)
            except InfeasibleLabels as err:
                raise InvalidSpatialCrystal(str(err),err.witness) from err
        self.last_dx, self.last_labels, self.last_info = dx.copy(), labels, orientation
        return np.r_[self.bath.energy(dx),crystal,harmonic.gauge_energy(labels)].astype(np.longdouble)


class HarmonicShapeDriver(Driver):
    """Keep the same melt gauge for active-label trials, signed shapes and accepted shapes."""
    def parts_at(self, t, spatial):
        return ReducedShapeWork(self,t,spatial)

    def step(self, force, curvature, parts, energy, t):
        previous = self.x.copy()
        bound = super().step(force,curvature,parts,energy,t)
        if parts.last_dx is None or not np.array_equal(self.x,previous+parts.last_dx):
            raise ArithmeticError('Accepted shape and last feasible slaved-label trial differ')
        self.labels = parts.last_labels
        return bound
