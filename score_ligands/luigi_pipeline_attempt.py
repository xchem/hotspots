import luigi
import os
from glob import glob
from os.path import join
from DiamondRunner import DiamondRunner
from DiamondScorer import DiamondScorer


class LDiamondRunner(luigi.Task):

    input_dict = luigi.parameter.DictParameter()

    def run(self):
        ldr = DiamondRunner(stem_dir=self.input_dict["stem"],
                            prot_name=self.input_dict["protein_name"],
                            xstal_id=self.input_dict["event_id"],
                            keyword=self.input_dict["keyword"],
                            chains=self.input_dict["chains"])
        ldr.run_hotspot_calculation(self.input_dict["number_rotations"])
        #ldr.extract_ligands()

    def output(self):
        path_list = glob(join(self.input_dict["stem"],
                             "hotspot_results",
                             "{}_{}_{}_*_results".format(self.input_dict["protein_name"], self.input_dict["event_id"], "".join(self.input_dict["chains"])),
                              "out.zip"))
        tar_list = [luigi.LocalTarget(p) for p in path_list if self.input_dict["keyword"] in p]
        return tar_list



class ParalleliseLDiamondRunner(luigi.WrapperTask):
    stem = luigi.parameter.Parameter()
    protein_name = luigi.parameter.Parameter()
    keyword = luigi.parameter.Parameter()
    chains = luigi.parameter.ListParameter()
    number_rotations = luigi.parameter.IntParameter(default=100000)

    def _get_inputs(self):
        paths = glob(join(self.stem, "x*"))
        event_ids = [p.split("x")[1] for p in paths]
        return [{"stem": self.stem,
                 "protein_name": self.protein_name,
                 "event_id": ev_id,
                 "keyword": self.keyword,
                 "chains":self.chains,
                 "number_rotations": self.number_rotations} for ev_id in event_ids]

    def requires(self):
        for k in self._get_inputs():
            yield LDiamondRunner(k)




if __name__ == "__main__":
    prot_name = "NUDT5"
    #stem = join(os.getcwd(), prot_name)
    stem= "/dls/science/groups/i04-1/software/mihaela/Data/NUDT5A/NUDT5A_fragment_hits"
    lchains=["A", "B"]
    keyword="bound"
    #nrot = 100000
    luigi.build([ParalleliseLDiamondRunner(stem, prot_name, keyword, lchains)], local_scheduler=True, workers=7)