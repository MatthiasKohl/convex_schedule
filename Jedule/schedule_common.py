
# SCHEDULING INDEPENDENT MOLDABLE TASKS ON MULTI-CORES WITH GPUS
# Copyright (C) 2014 Sascha Hunold <sascha@hunoldscience.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from jinja2 import Environment, FileSystemLoader
import os

class MoldSchedule:

    def __init__(self, nb_cpu):
        self.task_rects = []
        self.nb_cpu = nb_cpu
        self.metainfo = {}

    def get_nb_cpu(self):
        return self.nb_cpu

    def add_task_rect(self, task_rect):
        self.task_rects.append(task_rect)

    def get_task_rects(self):
        return self.task_rects

    def get_makespan(self):
        makespan = 0.0
        for task_rect in self.task_rects:
            makespan = max(makespan, task_rect.get_end_time())
        return makespan

    def set_metainfo(self, key, value):
        self.metainfo[key] = value

    def get_metainfo(self, key):
        return self.metainfo[key]

    def get_jedule_output(self):
        template_dir = os.path.dirname(os.path.abspath(__file__))
        j2_env = Environment(
            loader=FileSystemLoader(template_dir),
            trim_blocks=True,
            lstrip_blocks=True
        )
        template = j2_env.get_template('jinja.jed')
        return template.stream(vars(self))


class TaskRect:

    def __init__(self, task_id, node_type, device_id=0):
        self.task_id = task_id
        self.node_type = node_type
        self.device_id = device_id
        self.resources = [] # list of contiguous blocs: (start_id, nb_procs)
        self.start_time = None
        self.end_time = None

    def get_task_id(self):
        return self.task_id

    def get_device_id(self):
        return self.device_id

    def get_nbp(self):
        # TODO: could be cached or computed when pus are added
        return sum(r[1] for r in self.resources)

    def get_type(self):
        return self.node_type

    def get_start_time(self):
        return self.start_time

    def get_end_time(self):
        return self.end_time

    def set_procs(self, resources):
        self.resources = resources

    def set_times(self, start_time, end_time):
        self.start_time = start_time
        self.end_time = end_time

    def print_rect(self):
        print("*********************************************")
        print("task id" +  str(self.task_id) + " device:" + str(self.device_id))
        print("resources (start_id, nb_pus)" + str(self.resources))
        print("start:" + str(self.start_time))
        print("end  :" + str(self.end_time))
