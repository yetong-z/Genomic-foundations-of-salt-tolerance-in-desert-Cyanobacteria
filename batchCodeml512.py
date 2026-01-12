
class Factory(object):

    @classmethod
    def creat_dir(self, folder):  
        if not os.path.isdir(folder):
            # os.mkdirs(folder)
            os.makedirs(folder, exist_ok=True)

    @classmethod
    def moveFiles(self, files, targetPath):
        
        for k in files:
            shutil.copy(k, targetPath)

    @classmethod
    def remove_dir(self, path): 
        filelist = os.listdir(path)
        for f in filelist:
            filepath = os.path.join(path, f)
            if os.path.isfile(filepath):
                os.remove(filepath)
                #print(filepath+" removed!")
            elif os.path.isdir(filepath):
                shutil.rmtree(filepath, True)

    @classmethod
    def file_not_empty(self, file):
        return (os.path.exists(file) and (os.path.getsize(file) > 0))


class CodemlPack(object):

    def __init__(self, **kargs):
        self.dictArgsOriginal = copy.deepcopy(kargs)
        self.dictArgs = kargs
        # self.dictArgsOriginal = copy.deepcopy(self.dictArgs)
        self.implement_revision = self.dictArgs["rev"]
        self.null = self.dictArgs["null"]
        self.codeml = dict_codeml_ctl["codeml"]
        self.InputTree = dict_codeml_ctl["InputTree"]
        self.InputFile = dict_codeml_ctl["InputFile"]
        self.workDir = dict_codeml_ctl["workDir"]

    def exec_(self):
        if self.implement_revision and self.null:
            self.revision_null("revision")
            self.parseMLC()
            self.add_node_label()
            self.revision_null("null")
            self.run()
            os.remove(self.tree_revision_path)
           
            # self.calculate_p_value()
#         self.runM0()
#         self.runBM()
       
        elif self.null:
            self.revision_null("null")
            self.run()
            
            # self.calculate_p_value()
        else:
            
            self.runM0()

    def revision_null_sub(self, flag="M0"):
        if flag == "M0":
            self.dictArgs["model"] = "0"
            self.dictArgs["NSsites"] = "0"
            self.dictArgs["fix_kappa"] = "0"
            self.dictArgs["kappa"] = "2"
            self.dictArgs["fix_omega"] = "0"
            self.dictArgs["omega"] = "2"
            self.dictArgs["fix_blength"] = "0"
            self.dictArgs["outfile"] = "MO.mlc"
        elif flag == "MAnull":
            self.dictArgs["model"] = "2"
            self.dictArgs["NSsites"] = "2"
            self.dictArgs["fix_kappa"] = "0"
            self.dictArgs["kappa"] = "2"
            self.dictArgs["fix_omega"] = "1"
            self.dictArgs["omega"] = "1"
            self.dictArgs["fix_blength"] = "0"
            self.dictArgs["outfile"] = "MAnull.mlc"
        elif flag == "M1a":
            self.dictArgs["model"] = "0"
            self.dictArgs["NSsites"] = "1"
            self.dictArgs["fix_kappa"] = "0"
            self.dictArgs["kappa"] = "2"
            self.dictArgs["fix_omega"] = "0"
            self.dictArgs["omega"] = "2"
            self.dictArgs["fix_blength"] = "0"
            self.dictArgs["outfile"] = "M1a.mlc"

    def revision_null(self, mode="revision"):
       
        if mode=="revision":
            self.revision_null_sub(self.implement_revision)
            self.workPath = self.workDir + os.sep + \
                            "revision"
            treePath = self.InputTree
        elif mode=="null":
        
            self.dictArgs = copy.deepcopy(self.dictArgsOriginal)
            self.revision_null_sub(self.null)
            self.workPath = self.workDir + \
                            os.sep + "null"
            if self.implement_revision:
                self.dictArgs["fix_kappa"] = "1"
                self.dictArgs["kappa"] = self.kappa
                self.dictArgs["fix_blength"] = "2"
                treePath = self.tree_revision_path if self.implement_revision else self.InputTree
                self.dictArgs["treefile"] = os.path.basename(treePath)
        self.ctl = self.refreshCTL()
        Factory.creat_dir(self.workPath)
        
        self.saveCTL(self.ctl, self.workPath)
        
        Factory.moveFiles(
            [self.InputFile, treePath, self.codeml], self.workPath)
        self.null_result = self.workPath + os.sep + self.dictArgs["outfile"]
        # print("###%s revision###" % self.implement_revision)
        self.runCodeml(self.workPath)
        # print("###%s revision finished###" % self.implement_revision)

    def parseMLC(self):
        
        with open(self.null_result) as f:
            mlc = f.read()
        rgx_kappa = re.compile(r"kappa \(ts\/tv\) \= +?(.+)\n\n")
        rgx_tree = re.compile(
            r"\n(.+?)\n\nDetailed output identifying parameters")
        try:
            self.kappa = rgx_kappa.search(mlc).group(1)
            self.tree = rgx_tree.search(mlc).group(1)
        except:
            print(self.null_result, "parse MLC failed")

    def add_node_label(self):
        from ete3 import Tree
        
        t1 = Tree(self.InputTree, format=1)
        t2 = Tree(self.tree, format=1)
        
        list_leaf_label = []
        
        dict_ = {} # {"a": ["#1", ["a]], "b&%&c": ["$1", ["b", "c"]]}
        for node in t1.traverse("postorder"):
            # Do some analysis on node
            if node.name in ["$1", "#1"]:
                list_leaves = [i.name for i in node.get_leaves()]
                dict_["&%&".join(list_leaves)] = [node.name, list_leaves]
            if node.is_leaf() and ("#1" in node.name or "$1" in node.name):
                list_leaf_label.append(node.name)
        treeBase = os.path.splitext(self.dictArgs["treefile"])[0]
        self.tree_revision_path = self.workDir + \
            os.sep + "%s_revision.nwk" % treeBase
        if dict_:
            for i in dict_:
                mark, list_leaves = dict_[i]
                ancestor = t2.get_common_ancestor(list_leaves)
                ancestor.name = mark
        if list_leaf_label != []:
            for j in list_leaf_label:
                for node in t2.traverse("postorder"):
                    if node.is_leaf() and node.name in j:
                        node.name = j
        t2.write(
            format=1, outfile=self.tree_revision_path)

    def run(self):
        treePath = self.tree_revision_path if self.implement_revision else self.InputTree
        self.dictArgs = copy.deepcopy(self.dictArgsOriginal)
        if self.implement_revision:
            self.dictArgs["fix_kappa"] = "1"
            self.dictArgs["kappa"] = self.kappa
            self.dictArgs["fix_blength"] = "2"
       
        self.dictArgs["treefile"] = os.path.basename(treePath)
        self.ctl = self.refreshCTL()
        
        self.altn_workPath = self.workDir + os.sep + "alternative"
        Factory.creat_dir(self.altn_workPath)
        
        self.saveCTL(self.ctl, self.altn_workPath)
        
        Factory.moveFiles(
            [self.codeml, treePath, self.InputFile], self.altn_workPath)
        print("###alternative model###")
        self.runCodeml(self.altn_workPath)
        print("###alternative model finished###")

    def runM0(self):
        treePath = self.tree_revision_path if self.implement_revision else self.InputTree
        if self.implement_revision:
            self.dictArgs["fix_kappa"] = "1"
            self.dictArgs["kappa"] = self.kappa
            self.dictArgs["fix_blength"] = "2"
        
        self.dictArgs["treefile"] = os.path.basename(treePath)
        self.dictArgs["model"] = "0"
        self.dictArgs["NSsites"] = "0"
        
        self.dictArgs["outfile"] = self.dictArgsOriginal["outfile"]
        self.ctl = self.refreshCTL()
        self.M0_workPath = self.workDir + os.sep + "M0"
        Factory.creat_dir(self.M0_workPath)
        
        self.saveCTL(self.ctl, self.M0_workPath)
        
        Factory.moveFiles(
            [self.codeml, treePath, self.InputFile], self.M0_workPath)
        print("###M0 model###")
        self.runCodeml(self.M0_workPath)
        print("###M0 model finished###")

    def runBM(self):
        treePath = self.tree_revision_path if self.implement_revision else self.InputTree
        self.dictArgs = self.dictArgsOriginal
        if self.implement_revision:
            self.dictArgs["fix_kappa"] = "1"
            self.dictArgs["kappa"] = self.kappa
            self.dictArgs["fix_blength"] = "2"
        
        self.dictArgs["treefile"] = os.path.basename(treePath)
        self.ctl = self.refreshCTL()
        
        self.BM_workPath = self.workDir + os.sep + "BM"
        Factory.creat_dir(self.BM_workPath)
        
        self.saveCTL(self.ctl, self.BM_workPath)
        
        Factory.moveFiles(
            [self.codeml, treePath, self.InputFile], self.BM_workPath)
        print("###BM model###")
        self.runCodeml(self.BM_workPath)
        print("###BM model finished###")

    def calculate_p_value(self, alter_mlc, null_mlc):
        
        with open(null_mlc) as f: # self.null_result
            rev_mlc = f.read()
        
        with open(alter_mlc) as f1: # self.altn_workPath + os.sep + mlc_file_name
            altn_mlc = f1.read()
        rgx_p_value = re.compile(r"lnL\(ntime:.*np: *?(\d+?)\): *?([-\d\.]+)")
       
        null_np, null_lnl = rgx_p_value.search(rev_mlc).groups()
        altn_np, altn_lnl = rgx_p_value.search(altn_mlc).groups()
#         altn_b_omega, self.altn_f_omega = rgx_omega.search(altn_mlc).groups()
        differ_lnl = 2 * abs(float(null_lnl) - float(altn_lnl))
        differ_np = abs(int(null_np) - int(altn_np))
        # chisqprob is deprecated! stats.chisqprob is deprecated in scipy 0.17.0; use stats.distributions.chi2.sf instead.
        self.P_value = str(chi2.sf(differ_lnl, differ_np))
        return self.P_value

    def runCodeml(self, workPath):
        
        os.chdir(workPath)
        self.command = "codeml.exe" if "Windows" in platform.platform(
        ) else "codeml"
        popen = subprocess.Popen(
            self.command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)  
        stdout = self.run_command(
                        self.command, popen)
        with open("log.txt", "w") as f:
            f.write(stdout)

    def run_command(self, commands, popen):  
        print(commands)
        stdout = []
        print("###working on tree: %s, file: %s###" %
              (self.dictArgs["treefile"], self.dictArgs["seqfile"]))
        while True:
            try:
                out_line = popen.stdout.readline().decode("utf-8")
            except UnicodeDecodeError:
                out_line = popen.stdout.readline().decode("gbk")
            if out_line == "" and popen.poll() is not None:
                break
            # print(out_line)
            stdout.append(out_line)
        return "".join(stdout)

    def refreshCTL(self):
        ctl = '''      seqfile = {self.dictArgs[seqfile]}
     treefile = {self.dictArgs[treefile]}
      outfile = {self.dictArgs[outfile]}

        noisy = {self.dictArgs[noisy]}   * 0,1,2,3,9: how much rubbish on the screen
      verbose = {self.dictArgs[verbose]}   * 1: detailed output, 0: concise output
      runmode = {self.dictArgs[runmode]}   * 0: user tree;  1: semi-automatic;  2: automatic
                    * 3: StepwiseAddition; (4,5):PerturbationNNI 

      seqtype = {self.dictArgs[seqtype]}   * 1:codons; 2:AAs; 3:codons-->AAs
    CodonFreq = {self.dictArgs[CodonFreq]}   * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
        clock = {self.dictArgs[clock]}   * 0: no clock, unrooted tree, 1: clock, rooted tree
       aaDist = {self.dictArgs[aaDist]}
        model = {self.dictArgs[model]}
                    * models for codons:
                        * 0:one, 1:b, 2:2 or more dN/dS ratios for branches

      NSsites = {self.dictArgs[NSsites]}   * 0:one w;1:neutral;2:selection;3:discrete;4:freqs;
                    * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
                    * 10:beta&gamma+1;11:beta&normal>1;12:0&2normal>1;
                    * 13:3normal>0
        icode = {self.dictArgs[icode]}   * 0:standard genetic code; 1:mammalian mt; 2-10:see below
        Mgene = {self.dictArgs[Mgene]}

    fix_kappa = {self.dictArgs[fix_kappa]}   * 1: kappa fixed, 0: kappa to be estimated
        kappa = {self.dictArgs[kappa]}   * initial or fixed kappa
    fix_omega = {self.dictArgs[fix_omega]}   * 1: omega or omega_1 fixed, 0: estimate 
        omega = {self.dictArgs[omega]}   * initial or fixed omega, for codons or codon-transltd AAs

    fix_alpha = {self.dictArgs[fix_alpha]}   * 0: estimate gamma shape parameter; 1: fix it at alpha
        alpha = {self.dictArgs[alpha]}  * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = {self.dictArgs[Malpha]}   * different alphas for genes
        ncatG = {self.dictArgs[ncatG]}   * # of categories in the dG or AdG models of rates

        getSE = {self.dictArgs[getSE]}   * 0: don't want them, 1: want S.E.s of estimates
 RateAncestor = {self.dictArgs[RateAncestor]}   * (1/0): rates (alpha>0) or ancestral states (alpha=0)

  fix_blength = {self.dictArgs[fix_blength]}  * 0: ignore, -1: random, 1: initial, 2: fixed
       method = {self.dictArgs[method]}   * 0: simultaneous; 1: one branch at a time
   Small_Diff = {self.dictArgs[Small_Diff]}
    cleandata = {self.dictArgs[cleandata]}

* Specifications for duplicating results for the small data set in table 1
* of Yang (1998 MBE 15:568-573).
* see the tree file lysozyme.trees for specification of node (branch) labels
'''.format(self=self)
        return ctl

    def saveCTL(self, content, path):
        with open(path + os.sep + "codeml.ctl", "w") as f:
            f.write(self.ctl)

if __name__ == "__main__":
    import re
    import os
    import sys
    import argparse
    import glob
    import shutil
    import subprocess
    import platform
    import copy
    # from scipy.stats import chisqprob #chisqprob is deprecated! stats.chisqprob is deprecated in scipy 0.17.0; use stats.distributions.chi2.sf instead.
    from collections import OrderedDict
    import multiprocessing
    Pool = multiprocessing.Pool
    from scipy.stats.distributions import chi2

    
    scriptPath = os.path.dirname(os.path.realpath(__file__))

    
    listFiles = glob.glob(scriptPath + os.sep + "InputFiles" + os.sep + "*")
    listTrees = glob.glob(scriptPath + os.sep + "InputTrees" + os.sep + "*")

    def parameter():
        parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            prog='batchCodeml.py',
            description='',
            epilog=r'''python batchCodeml.py /home/zhangdong/PAML/paml4.9e/bin/codeml -icode 6 -model 2 -NSsites 0 -omega 1 -rev M0
              ''')
        parser.add_argument(
            "codeml", help='excutable codeml path')
        parser.add_argument(
            '-f', dest='files', default=listFiles, nargs='*', help='input files')
        parser.add_argument(
            '-t', dest='trees', default=listTrees, nargs='*', help='input trees')
        parser.add_argument(
            '-noisy', dest='noisy', default="9")
        parser.add_argument(
            '-verbose', dest='verbose', default="0")
        parser.add_argument(
            '-runmode', dest='runmode', default="0")
        parser.add_argument(
            '-seqtype', dest='seqtype', default="1", help='1:codons; 2:AAs; 3:codons-->AAs')
        parser.add_argument(
            '-CodonFreq', dest='CodonFreq', default="2", help='0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table')
        parser.add_argument(
            '-clock', dest='clock', default="0", help='0: no clock, unrooted tree, 1: clock, rooted tree')
        parser.add_argument(
            '-aaDist', dest='aaDist', default="0")
        parser.add_argument(
            '-model', dest='model', default="0", help='0:one, 1:b, 2:2 or more dN/dS ratios for branches')
        parser.add_argument(
            '-NSsites', dest='NSsites', default="0", help='0:one w;1:neutral;2:selection;3:discrete;4:freqs;5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;10:beta&gamma+1;11:beta&normal>1;12:0&2normal>1;13:3normal>0')
        parser.add_argument(
            '-icode', dest='icode', default="0")
        parser.add_argument(
            '-Mgene', dest='Mgene', default="0")
        parser.add_argument(
            '-fix_kappa', dest='fix_kappa', default="0")
        parser.add_argument(
            '-kappa', dest='kappa', default="0")
        parser.add_argument(
            '-fix_omega', dest='fix_omega', default="0")
        parser.add_argument(
            '-omega', dest='omega', default="2")
        parser.add_argument(
            '-alpha', dest='alpha', default=".0", help='when activate estimate just type "estimate"')
        parser.add_argument(
            '-Malpha', dest='Malpha', default="0")
        parser.add_argument(
            '-ncatG', dest='ncatG', default="3")
        parser.add_argument(
            '-getSE', dest='getSE', default="0")
        parser.add_argument(
            '-RateAncestor', dest='RateAncestor', default="0")
        parser.add_argument(
            '-fix_blength', dest='fix_blength', default="0")
        parser.add_argument(
            '-method', dest='method', default="0")
        parser.add_argument(
            '-Small_Diff', dest='Small_Diff', default=".45e-6")
        parser.add_argument(
            '-cleandata', dest='cleandata', default="0")
        parser.add_argument('-rev', dest='rev', help='implement revision',
                            default=None, choices=["M0", "MAnull", "M1a", None])
        parser.add_argument('-null', dest='null', help='null',
                            default=None, choices=["M0", "MAnull", "M1a", None])
        parser.add_argument('-cpu', dest='threads', help='the threads for parallel run',
            default=10, type=int)
        myargs = parser.parse_args(sys.argv[1:])
        return myargs
    myargs = parameter()
    dict_codeml_ctl = OrderedDict()
    dict_codeml_ctl["noisy"] = myargs.noisy
    dict_codeml_ctl["verbose"] = myargs.verbose
    dict_codeml_ctl["runmode"] = myargs.runmode
    dict_codeml_ctl["seqtype"] = myargs.seqtype
    dict_codeml_ctl["CodonFreq"] = myargs.CodonFreq
    dict_codeml_ctl["clock"] = myargs.clock
    dict_codeml_ctl["aaDist"] = myargs.aaDist
    dict_codeml_ctl["model"] = myargs.model
    dict_codeml_ctl["NSsites"] = myargs.NSsites
    dict_codeml_ctl["icode"] = myargs.icode
    dict_codeml_ctl["Mgene"] = myargs.Mgene
    dict_codeml_ctl["fix_kappa"] = myargs.fix_kappa
    dict_codeml_ctl["kappa"] = myargs.kappa
    dict_codeml_ctl["fix_omega"] = myargs.fix_omega
    dict_codeml_ctl["omega"] = myargs.omega
    dict_codeml_ctl["fix_alpha"] = "0" if myargs.alpha == "estimate" else "1"
    dict_codeml_ctl[
        "alpha"] = myargs.alpha if myargs.alpha != "estimate" else "0"
    dict_codeml_ctl["Malpha"] = myargs.Malpha
    dict_codeml_ctl["ncatG"] = myargs.ncatG
    dict_codeml_ctl["getSE"] = myargs.getSE
    dict_codeml_ctl["RateAncestor"] = myargs.RateAncestor
    dict_codeml_ctl["fix_blength"] = myargs.fix_blength
    dict_codeml_ctl["method"] = myargs.method
    dict_codeml_ctl["Small_Diff"] = myargs.Small_Diff
    dict_codeml_ctl["cleandata"] = myargs.cleandata

    outputDir = scriptPath + os.sep + "Results"
    Factory.creat_dir(outputDir)
    # Factory.remove_dir(outputDir)
    
    codeml_absPath = os.path.abspath(myargs.codeml)
    sum = "File,Tree,P-value\n"
    print(" ".join(sys.argv))
    print(f"Totally {len(myargs.trees)} trees and {len(myargs.files)} files")
    def run_(i, j):
        
        tree_base_prefix = os.path.splitext(os.path.basename(i))[0]
        workDir = scriptPath + os.sep + "Results" + os.sep + tree_base_prefix
        Factory.creat_dir(workDir)
        file_base_prefix = os.path.splitext(os.path.basename(j))[0]
        dict_codeml_ctl["seqfile"] = os.path.basename(j)
        dict_codeml_ctl["treefile"] = os.path.basename(i)
        dict_codeml_ctl["outfile"] = file_base_prefix + ".mlc"
        dict_codeml_ctl["codeml"] = codeml_absPath
        dict_codeml_ctl["InputTree"] = i
        dict_codeml_ctl["InputFile"] = j
        dict_codeml_ctl["rev"] = myargs.rev
        dict_codeml_ctl["null"] = myargs.null
        
        workDir_each = workDir + os.sep + file_base_prefix
        Factory.creat_dir(workDir_each)
        dict_codeml_ctl["workDir"] = workDir_each
        alter_mlc = f"{workDir_each}/alternative/{dict_codeml_ctl['outfile']}"
        # print(alter_mlc)
        if Factory.file_not_empty(alter_mlc):
            print("###The result for tree: %s, file: %s already generated!###" %
                  (tree_base_prefix, file_base_prefix))
            codemlpack = CodemlPack(**dict_codeml_ctl)
            codemlpack.revision_null_sub(flag=dict_codeml_ctl["null"])
            null_mlc = f"{workDir_each}/null/{codemlpack.dictArgs['outfile']}"
            p = codemlpack.calculate_p_value(alter_mlc, null_mlc)
        else:
            print("###working on tree: %s, file: %s###" %
                  (tree_base_prefix, file_base_prefix))
            codemlpack = CodemlPack(**dict_codeml_ctl)
            codemlpack.exec_()
            p = codemlpack.calculate_p_value(codemlpack.altn_workPath + os.sep + codemlpack.dictArgs["outfile"],
                                              codemlpack.null_result)
        # print(file_base_prefix, tree_base_prefix, p)
        return ",".join([file_base_prefix, tree_base_prefix, p]) + "\n"
    
    threads = myargs.threads
    threads = threads if len(myargs.files) > threads else len(myargs.files)
    pool = Pool(processes=threads)
    try:
        r = [pool.apply_async(run_, (i, j)) for i in myargs.trees for j in myargs.files]
        pool.close()
        for item in r:
            item.wait(timeout=9999999)  # Without a timeout, you can't interrupt this.
            sum += item.get()
    except KeyboardInterrupt:
        pool.terminate()
    finally:
        pool.join()
    #
    # for j in myargs.files:
    #     file_base_prefix = os.path.splitext(os.path.basename(j))[0]
    #     dict_codeml_ctl["seqfile"] = os.path.basename(j)
    #     dict_codeml_ctl["treefile"] = os.path.basename(i)
    #     dict_codeml_ctl["outfile"] = file_base_prefix + ".mlc"
    #     dict_codeml_ctl["codeml"] = codeml_absPath
    #     dict_codeml_ctl["InputTree"] = i
    #     dict_codeml_ctl["InputFile"] = j
    #     dict_codeml_ctl["rev"] = myargs.rev
    #     dict_codeml_ctl["null"] = myargs.null
    #     workDir_each = workDir + os.sep + file_base_prefix
    #     Factory.creat_dir(workDir_each)
    #     dict_codeml_ctl["workDir"] = workDir_each
    #     print("###working on tree: %s, file: %s###" %
    #           (tree_base_prefix, file_base_prefix))
    #     codemlpack = CodemlPack(**dict_codeml_ctl)
    #     sum += ",".join([file_base_prefix, tree_base_prefix,
    #                      codemlpack.P_value]) + "\n"
    with open(outputDir + os.sep + "sum.csv", "w") as f:
        f.write(sum)

    print("Done!")

