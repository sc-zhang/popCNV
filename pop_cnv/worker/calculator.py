from pathos.multiprocessing import Pool
from pop_cnv.io.loader import FastaLoader
from pop_cnv.io.message import Message
from pop_cnv.worker.generator import BIN
from pop_cnv.worker.runner import Runner
from os import listdir, path, chdir
from numpy import median, isnan, isinf, bincount
from scipy.stats import f_oneway
import pysam
import gzip


class GC:
    def __init__(self):
        self.gc_db = {}
        self.__bins = {}
        self.__fa_db = {}

    @staticmethod
    def __sub_stat(bins, seq, chrn):
        gc_list = []
        for sp, ep in bins:
            gc_cnt = 0
            for i in range(sp - 1, ep):
                if seq[i].lower() == 'g' or seq[i].lower() == 'c':
                    gc_cnt += 1
            gc_list.append([sp, ep, int(gc_cnt * 100.0 / (ep - sp + 1))])
        return gc_list, chrn

    def stat(self, genome_file, window_size, threads):
        self.__fa_db = FastaLoader.load(genome_file)
        self.__bins = BIN.generate(window_size, self.__fa_db)

        pool = Pool(processes=threads)
        res = []
        for chrn in self.__fa_db:
            r = pool.apply_async(self.__sub_stat, args=(self.__bins, self.__fa_db[chrn], chrn, ))
            res.append(r)
        pool.close()
        pool.join()

        for r in res:
            gc_list, chrn = r.get()
            for sp, ep, gc_prop in gc_list:
                self.gc_db[tuple([chrn, sp, ep])] = gc_prop


class SeqDepth:
    def __init__(self):
        self.depth_db = {}

    @staticmethod
    def __sub_calc(genome_size, bam_file, smp):
        read_size = 0
        with pysam.AlignmentFile(bam_file, 'rb') as fin:
            for line in fin:
                ql = line.query_length
                if ql != -1:
                    read_size += ql
        return read_size * 1.0 / genome_size, smp

    def calc(self, genome_file, bam_path, threads):
        genome_size = 0
        fa_db = FastaLoader.load(genome_file)
        for _ in fa_db:
            genome_size += len(fa_db[_])

        pool = Pool(processes=threads)
        res = []
        for bam in listdir(bam_path):
            if bam.split('.')[-1].lower() != 'bam':
                continue
            smp = bam.split('.')[0]
            r = pool.apply_async(self.__sub_calc, args=(genome_size, path.join(bam_path, bam), smp,))
            res.append(r)
        pool.close()
        pool.join()

        for r in res:
            depth, smp = r.get()
            self.depth_db[smp] = depth


class BamDepth:
    def __init__(self):
        pass

    @staticmethod
    def __sub_run(sample, bam_file, bin_size, out_dir):
        runner = Runner()
        chdir(out_dir)
        cmd = "mosdepth -b %d %s %s" % (bin_size, sample, bam_file)
        runner.set_cmd(cmd)
        runner.run()

    def run(self, bam_path, bin_size, out_dir, threads):
        pool = Pool(processes=threads)
        for bam in listdir(bam_path):
            if bam.split('.')[-1].lower() != 'bam':
                continue
            smp = bam.split('.')[0]
            bam_file = path.join(bam_path, bam)
            pool.apply_async(self.__sub_run, args=(smp, bam_file, bin_size, out_dir,))
        pool.close()
        pool.join()


class Norm:
    def __init__(self):
        self.norm_db = {}

    @staticmethod
    def __sub_norm(mos_file, smp, s_ratio):
        rd_db = {}
        with gzip.open(mos_file, 'rt') as fin:
            for line in fin:
                data = line.strip().split()
                chrn = data[0]
                sp = int(data[1])+1
                ep = int(data[2])
                rd = float(data[-1])*1./s_ratio
                rd_db[tuple([chrn, sp, ep])] = rd
        return rd_db, smp

    def norm(self, mos_path, sample_depth_db, threads):
        pool = Pool(processes=threads)
        res = []
        for fn in listdir(mos_path):
            if not fn.endswith(".regions.bed.gz"):
                continue
            smp = fn.split('.')[0]
            if smp not in sample_depth_db:
                continue
            mos_file = path.join(mos_path, fn)
            r = pool.apply_async(self.__sub_norm, args=(mos_file, smp, sample_depth_db[smp],))
            res.append(r)
        pool.close()
        pool.join()

        for r in res:
            rd_db, smp = r.get()
            self.norm_db[smp] = rd_db


class CN:
    def __init__(self):
        self.cn_db = {}

    @staticmethod
    def __sub_convert(gc_db, sub_rd_db, smp):
        conv_gc_db = {}
        for _ in gc_db:
            gc = gc_db[_]
            rd = sub_rd_db[_]
            if gc not in conv_gc_db:
                conv_gc_db[gc] = []
            conv_gc_db[gc].append(rd)

        conv_cn_db = {}
        for gc in conv_gc_db:
            cn = median(gc_db[gc])
            conv_cn_db[gc] = cn
        cn_db = {}
        for _ in gc_db:
            cn_db[_] = sub_rd_db[_]*1./conv_cn_db[gc_db[_]]
        return cn_db, smp

    def convert(self, gc_db, rd_db, threads):
        pool = Pool(processes=threads)
        res = []
        for smp in rd_db:
            r = pool.apply_async(self.__sub_convert, args=(gc_db, rd_db[smp], smp, ))
            res.append(r)
        pool.close()
        pool.join()

        for r in res:
            cn_db, smp = r.get()
            self.cn_db[smp] = cn_db


class GeneCN:
    def __init__(self):
        self.gene_cn = {}

    @staticmethod
    def search_pos(pos_list, pos):
        s = 0
        e = len(pos_list) - 1
        while s <= e:
            mid = int((s + e) / 2)
            if pos_list[mid][0] > pos:
                e = mid - 1
            elif pos_list[mid][0] < pos:
                s = mid + 1
            else:
                return mid
        if pos_list[e][1] >= pos:
            return e
        else:
            return -1

    def calc(self, cn_db, gene_bed):
        conv_cn_db = {}
        for smp in cn_db:
            conv_cn_db[smp] = {}
            for _ in cn_db[smp]:
                chrn, sp, ep = _
                if chrn not in conv_cn_db[smp]:
                    conv_cn_db[smp][chrn] = []
                conv_cn_db[smp][chrn].append([sp, ep, cn_db[smp][_]])
        gn_cn = {}
        for gn in gene_bed:
            gn_cn[gn] = {}
            for smp in conv_cn_db:
                chrn, sp, ep = gene_bed[gn]
                nsp = self.search_pos(conv_cn_db[smp][chrn], sp)
                nep = self.search_pos(conv_cn_db[smp][chrn], ep)
                if nsp == -1 or nep == -1:
                    continue
                cn_list = []
                for i in range(nsp, nep + 1):
                    cn_list.append(conv_cn_db[smp][chrn][i][2])
                try:
                    gn_cn[gn][smp] = median(cn_list)
                except Exception as e:
                    Message().warn(repr(e))
                    gn_cn[gn][smp] = float('nan')
        self.gene_cn = gn_cn


class RoundCN:
    def __init__(self):
        self.round_cn = {}

    def round(self, gn_cn):
        for gn in gn_cn:
            self.round_cn[gn] = {}
            for smp in gn_cn[gn]:
                self.round_cn[gn][smp] = int(round(gn_cn[gn][smp]))


class RFD:
    def __init__(self):
        self.rfd_db = {}

    def calc(self, round_cn, grp_db, wild_grp):
        for gn in round_cn:
            cn_db = {}
            pop_list = []
            for smp in round_cn[gn]:
                if isnan(round_cn[gn][smp]) or isinf(round_cn[gn][smp]):
                    continue
                grp = grp_db[smp]
                if grp not in cn_db:
                    cn_db[grp] = []
                cn_db[grp].append(round_cn[gn][smp])
                pop_list.append(round_cn[gn][smp])
            f_db = {}
            if len(pop_list) == 0:
                continue
            f_pop = max(bincount(pop_list))*1./len(pop_list)
            for grp in cn_db:
                if len(cn_db[grp]) == 0:
                    continue
                f_db[grp] = max(bincount(cn_db[grp]))*1./len(cn_db[grp])
            if wild_grp not in f_db:
                continue
            for grp in f_db:
                if grp == wild_grp:
                    continue
                _, p = f_oneway(cn_db[grp], cn_db[wild_grp])

                RFD = (f_db[grp]-f_db[wild_grp])*1./f_pop
                if grp not in self.rfd_db:
                    self.rfd_db[grp] = []
                self.rfd_db[grp].append([RFD, p, gn])
