import luigi
import os
from glob import glob
from os.path import join
from DiamondRunner import DiamondRunner
from DiamondScorer import DiamondScorer
import datetime


class LDiamondRunner(luigi.Task):
    """
    Luigi wrapper for running Hotspots on Diamond
    """

    input_dict = luigi.parameter.DictParameter()

    def run(self):
        ldr = DiamondRunner(stem_dir=self.input_dict["stem"],
                            prot_name=self.input_dict["protein_name"],
                            xstal_id=self.input_dict["xstal_id"],
                            keyword=self.input_dict["keyword"],
                            chains=self.input_dict["chains"],
                            target_save_dir=self.input_dict["tar_dir"])
        ldr.run_hotspot_calculation(self.input_dict["number_rotations"])
        # ldr.extract_ligands()

    def output(self):
        #print(self.input_dict)
        path_list = glob(join(self.input_dict["tar_dir"],
                              "{}_{}_{}_*_results".format(self.input_dict["protein_name"],
                                                          self.input_dict["xstal_id"],
                                                          "".join(self.input_dict["chains"])),
                              "out.zip"))
        tar_list = [luigi.LocalTarget(p) for p in path_list if self.input_dict["keyword"] in p]
        return tar_list


class LDiamondScorer(luigi.Task):
    """
    Luigi wrapper for hotspot  scoring of XChem fragment hits.
    """
    input_dict = luigi.parameter.DictParameter()

    def requires(self):
        return LDiamondRunner(self.input_dict)

    def run(self):
        lds = DiamondScorer(stem=self.input_dict["tar_dir"],
                            prot_name=self.input_dict["protein_name"],
                            x_id=self.input_dict["xstal_id"],
                            keyword=self.input_dict["keyword"],
                            chains=self.input_dict["chains"])
        lds.score_ligands()
        self._log_LDS(write_file=True)

    def _log_LDS(self, write_file=False):
        now = datetime.datetime.now()
        if self.input_dict["keyword"]:
            fname = join(self.input_dict["tar_dir"],
                         "{}_{}_{}_{}_results".format(self.input_dict["protein_name"],
                                                      self.input_dict["xstal_id"],
                                                      "".join(self.input_dict["chains"]),
                                                      self.input_dict["keyword"]),
                         "LDiamondScorer_complete.log")
        else:
            fname = join(self.input_dict["tar_dir"],
                         "{}_{}_{}_results".format(self.input_dict["protein_name"],
                                                      self.input_dict["xstal_id"],
                                                      "".join(self.input_dict["chains"])),
                         "LDiamondScorer_complete.log")

        if write_file:
            with open(fname, "w") as f:
                f.write(
                    "LuigiDiamondScorer log calculated at {}-{}-{} {}:{} \n".format(now.year, now.month,now.day, now.hour, now.minute))
                for key, val in self.input_dict.items():
                    f.write(str(key) + ": " + str(val) + "\n")
        else:
            return fname

    def output(self):

        return luigi.LocalTarget(self._log_LDS())


class PipelineTask(luigi.WrapperTask):
    """
    Luigi wrapper for scoring and calculating mini-pipeline
    """
    xstals_dir = luigi.parameter.Parameter(default=os.getcwd())
    hs_results_dir = luigi.parameter.Parameter(default=os.path.join(os.getcwd(), "hotspot_results"))
    prot_name = luigi.parameter.Parameter(default="myprotein")
    keyword = luigi.parameter.OptionalParameter(default="")
    chains = luigi.parameter.ListParameter(default=["A"])
    number_rotations = luigi.parameter.IntParameter(default=100000)

    def _get_inputs(self):
        paths = glob(join(self.xstals_dir, "x*"))
        xstal_ids = [p.split("x")[1] for p in paths]
        if isinstance(self.chains[0], tuple):
            return [{"stem": self.xstals_dir,
                     "protein_name": self.prot_name,
                     "xstal_id": x_id,
                     "keyword": self.keyword,
                     "chains": sub_chains,
                     "tar_dir": self.hs_results_dir,
                     "number_rotations": self.number_rotations} for x_id in xstal_ids for sub_chains in self.chains]

        else:

            return [{"stem": self.xstals_dir,
                     "protein_name": self.prot_name,
                     "xstal_id": x_id,
                     "keyword": self.keyword,
                     "chains": self.chains,
                     "tar_dir":self.hs_results_dir,
                     "number_rotations": self.number_rotations} for x_id in xstal_ids]

    def requires(self):
        for m in self._get_inputs():
            yield LDiamondScorer(m)

    def run(self):
        """
        
        :return: 
        """
        now = datetime.datetime.now()
        with open(join(self.xstals_dir, "Pipeline_complete.log"), "w") as f:
            f.write("Pipeline log calculated at {}-{}-{} {}:{} \n".format(now.year, now.month, now.day, now.hour, now.minute))
            for key, val in self.__dict__.items():
                f.write(str(key) + ": " + str(val) + "\n")

    # def output(self):
    #     return luigi.LocalTarget(join(self.xstals_dir, "Pipeline_complete.log"))


if __name__ == "__main__":
    pname = "PARP14"
    stem= "/home/jin76872/Desktop/Mih/Data/PARP14/aligned_hotspot_maps/aligned_structures/"
    res_dir = stem
    lchains = ["A"]
    #key = "bound"
    #nrot = 100000
    luigi.build([PipelineTask(xstals_dir=stem, prot_name=pname, keyword=None, chains=lchains, hs_results_dir=res_dir)], local_scheduler=True, workers=20)
