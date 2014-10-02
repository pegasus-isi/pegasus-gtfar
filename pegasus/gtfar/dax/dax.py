# Copyright 2007-2014 University Of Southern California
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

__author__ = 'Rajiv Mayani'

import os
import sys
import mmh3

from pegasus.gtfar.dax.AutoADAG import AutoADAG

from Pegasus.DAX3 import ADAG, Dependency, Job, File, Link, Executable, PFN, Transformation


class UNIXUtils(object):
    @staticmethod
    def cat(inputs, output, o_link=Link.OUTPUT, o_transfer=False, o_register=False):
        cat = Job(name='merge')

        # Outputs
        output = File(output)

        for input_file in inputs:
            # Inputs
            input_file = File(input_file)

            # Arguments
            cat.addArguments(input_file)

            # Uses
            cat.uses(input_file, link=Link.INPUT)

        cat.setStdout(output)
        cat.uses(output, link=o_link, transfer=o_transfer, register=o_register)

        return cat


class AnnotateMixin(object):
    def annotate(self):
        self._annotate(self._read_length)

        if self._is_map_filtered or self._is_trim_unmapped:
            for trim in self._trims:
                self._annotate(trim)

    def _annotate(self, read_length):
        seed = 'F%d' % self._compute_seed(read_length)

        self._annotate_gtf(read_length)

        self._features_index(read_length, seed=seed)
        self._chrs_index(read_length, seed=seed)
        self._splices_index(read_length, seed=seed)

        if self._clip_reads:
            # Seed is always F1
            self._genes_index(read_length, seed='F1')

    def _annotate_gtf(self, read_length):
        annotate_gtf = Job(name='annotate_gtf')
        annotate_gtf.invoke('all', self._state_update % 'Generating annotation FASTA files')

        prefix = self._get_index_hash(read_length, exclude_genome=True)

        # Inputs
        gtf = File(self._gtf)

        chromosomes = [str(i) for i in range(1, 23)]
        chromosomes.extend(['X', 'Y', 'R', 'M'])

        for i in chromosomes:
            chr_i = File('chr%s.fa' % i)

            # Uses
            annotate_gtf.uses(chr_i, link=Link.INPUT)

        # Outputs
        features = File('h%s/FEATURES.fa' % prefix)
        chrs = File('h%s/GENOME.fa' % prefix)
        splices = File('h%s/SPLICES.fa' % prefix)
        genes = File('h%s/GENE.fa' % prefix)

        # Arguments
        annotate_gtf.addArguments(gtf, '-c .', '-p h%s/' % prefix, '-l %d' % read_length)

        # Uses
        annotate_gtf.uses(gtf, link=Link.INPUT)
        annotate_gtf.uses(features, link=Link.OUTPUT, transfer=False, register=False)
        annotate_gtf.uses(chrs, link=Link.OUTPUT, transfer=False, register=False)
        annotate_gtf.uses(splices, link=Link.OUTPUT, transfer=False, register=False)
        annotate_gtf.uses(genes, link=Link.OUTPUT, transfer=False, register=False)

        self.adag.addJob(annotate_gtf)

    def _features_index(self, read_length, read_format='fastq', seed='F2'):
        return self._perm_index('FEATURES', read_length, read_format=read_format, seed=seed)

    def _chrs_index(self, read_length, read_format='fastq', seed='F2'):
        return self._perm_index('GENOME', read_length, read_format=read_format, seed=seed)

    def _splices_index(self, read_length, read_format='fastq', seed='F2'):
        return self._perm_index('SPLICES', read_length, read_format=read_format, seed=seed)

    def _genes_index(self, read_length, read_format='fastq', seed='F1'):
        return self._perm_index('GENE', read_length, read_format=read_format, seed=seed)

    def _perm_index(self, index_type, read_length, read_format='fastq', seed='F2'):
        perm_index = Job(name='perm')
        perm_index.invoke('all', self._state_update % 'Pre-computing %s index file' % index_type.capitalize())

        prefix = self._get_index_hash(read_length, exclude_genome=True)

        # Input files
        fa_input = File('h%s/%s.fa' % (prefix, index_type))

        # Output files
        hash_v = self._get_index_hash(read_length, seed)
        index = File('h%d_%s_%s_%s.index' % (hash_v, index_type, seed, read_length))

        # Arguments
        perm_index.addArguments(fa_input, '%d' % read_length, '--readFormat', read_format, '--seed', seed)
        perm_index.addArguments('-s', index)

        # Uses
        perm_index.uses(fa_input, link=Link.INPUT)
        # Save this file
        perm_index.uses(index, link=Link.OUTPUT, transfer=True, register=True)

        self.adag.addJob(perm_index)

        return perm_index


class FilterMixin(object):
    def option_filter(self):
        splits = self._splits
        suffix_len = max(2, len(str(splits)))
        rejects = []
        stats = []

        self._fastq_split(splits=splits, suffix_len=suffix_len)

        for i in range(splits):
            rejects.append('reads%d_reject.fastq' % i)
            stats.append('reads%d.stats' % i)
            self._pre_filter_fastq(i, suffix_len)

        # Merge rejects
        cat = UNIXUtils.cat(rejects, '%s.reject.fastq' % self._prefix, o_transfer=True)
        cat.invoke('all', self._state_update % 'Merge rejected reads file')
        self.adag.addJob(cat)

        self._merge_stats()

    def _fastq_split(self, splits=2, suffix_len=2):
        fastq_split = Job(name='fastq-split')
        fastq_split.invoke('all', self._state_update % 'Splitting input reads file into %d parts' % splits)

        # Inputs
        reads = File(self._reads)

        # Arguments
        fastq_split.addArguments(reads, '%d' % splits)

        # Uses
        fastq_split.uses(reads, link=Link.INPUT)

        for i in range(splits):
            split_i = File(('x%0' + str(suffix_len) + 'd') % i)

            # Outputs
            fastq_split.uses(split_i, link=Link.OUTPUT, transfer=False, register=False)

        self.adag.addJob(fastq_split)

    def _pre_filter_fastq(self, index, suffix_len):
        pre_filter = Job(name='pre_filter_fastq.py')
        pre_filter.invoke('all', self._state_update % 'Pre-filter reads part %d file' % (index + 1))
        prefix = 'reads%d' % index

        # Inputs
        reads = File(('x%0' + str(suffix_len) + 'd') % index)

        # Outputs
        full_fastq = File('%s_full.fastq' % prefix)
        reject = File('%s_reject.fastq' % prefix)
        stats = File('%s.stats' % prefix)

        # Arguments
        trims = ','.join([str(i) for i in self._trims])
        trims = '0' if trims == ',' else trims

        pre_filter.addArguments(reads, '-r', '%d' % self._read_length, '-t', '%s' % trims)
        pre_filter.addArguments('-p', prefix)

        # Uses
        pre_filter.uses(reads, link=Link.INPUT)

        for t in self._trims:
            fastq_t = File('%s_%d.fastq' % (prefix, t))
            pre_filter.uses(fastq_t, link=Link.OUTPUT, transfer=False, register=False)

        pre_filter.uses(full_fastq, link=Link.OUTPUT, transfer=False, register=False)
        pre_filter.uses(reject, link=Link.OUTPUT, transfer=False, register=False)
        pre_filter.uses(stats, link=Link.OUTPUT, transfer=False, register=False)

        self.adag.addJob(pre_filter)

    def _merge_stats(self):
        merge_stats = Job(name='merge-stats')
        merge_stats.invoke('all', self._state_update % 'Merging adaptor stats file')

        # Outputs
        adaptor_stats = File('%s.adaptor.stats' % self._prefix)

        # Arguments
        merge_stats.addArguments('reads*.stats', adaptor_stats)

        for i in range(self._splits):
            # Inputs
            stats_i = File('reads%d.stats' % i)

            # Uses
            merge_stats.uses(stats_i, link=Link.INPUT)

        # Outputs
        merge_stats.uses(adaptor_stats, link=Link.OUTPUT, transfer=True, register=False)

        self.adag.addJob(merge_stats)


class IterativeMapMixin(object):
    def iterative_map(self):
        if self._clip_reads:
            self._clip_and_parse_reads('reads%d_full_miss_FEATURES_miss_GENOME_miss_SPLICES.fastq', 'full')

        unmapped_reads = self._map_and_parse_reads('reads%d_full.fastq', 'full')

        for trim in self._trims:
            self._read_length = trim

            if self._is_trim_unmapped:
                unmapped_reads = self._map_and_parse_reads(unmapped_reads, 'trim%d' % trim)

            if self._is_map_filtered:
                self._map_and_parse_reads('reads%%d_%d.fastq' % trim, 'filter%d' % trim)

        cat = UNIXUtils.cat(self._vis_files, '%s.sam' % self._prefix, o_transfer=True)
        cat.invoke('all', self._state_update % 'Merge all vis files into a single SAM file')
        self.adag.addJob(cat)

        if self._is_trim_unmapped:
            cat = UNIXUtils.cat([2], '%s.unmapped.fastq')
        else:
            merged_unmapped = '%s.unmapped.fastq' % self._prefix
            cat = UNIXUtils.cat([unmapped_reads % i for i in self._range()], merged_unmapped, o_transfer=True)

        cat.invoke('all', self._state_update % 'Merge all unmapped reads')
        self.adag.addJob(cat)

    def _map_and_parse_reads(self, reads, tag):
        self._setup_perm_seeds()

        self._map_and_parse_reads_to_features(reads, tag)

        path, file_name, ext = GTFAR._get_filename_parts(reads)

        miss_features = '%s_miss_FEATURES%s' % (file_name, ext)
        self._map_and_parse_reads_to_genome(miss_features, tag)

        miss_genome = '%s_miss_GENOME%s' % os.path.splitext(miss_features)

        if self._splice:
            self._map_and_parse_reads_to_splices(miss_genome, tag)
            miss_splices = '%s_miss_SPLICES%s' % os.path.splitext(miss_genome)

            return miss_splices

        return miss_genome

    def _setup_perm_seeds(self):
        self._seed = self._mismatches

        if self._read_length > 64:
            self._seed = (self._mismatches + 1) / 2

    def _compute_seed(self, read_length):
        seed = self._mismatches

        if read_length > 64:
            seed = (self._mismatches + 1) / 2

        return seed

    def _map_and_parse_reads_to_features(self, reads, tag):
        self._perm('features', 'FEATURES', reads, tag)

        reads = os.path.basename(reads)
        for i in self._range():
            reads_i = os.path.splitext(reads)[0] % i
            mapping_file = 'FEATURES_B_%d_%d_%s.mapping' % (self._seed, self._mismatches, reads_i)
            self._parse_alignment(mapping_file, tag)

    def _map_and_parse_reads_to_genome(self, reads, tag):
        self._perm('chrs', 'GENOME', reads, tag, output_sam=True)

        reads = os.path.basename(reads)
        for i in self._range():
            reads_i = os.path.splitext(reads)[0] % i
            sam_file = 'GENOME_B_%d_%d_%s.sam' % (self._seed, self._mismatches, reads_i)
            self._parse_alignment(sam_file, tag)

    def _map_and_parse_reads_to_splices(self, reads, tag):
        self._perm('jxnCands', 'SPLICES', reads, tag)

        reads = os.path.basename(reads)
        for i in self._range():
            reads_i = os.path.splitext(reads)[0] % i
            mapping_file = 'SPLICES_B_%d_%d_%s.mapping' % (self._seed, self._mismatches, reads_i)
            self._parse_alignment(mapping_file, tag)

    def _perm(self, index_type, map_to, reads, tag, output_sam=False):
        perm = Job(name='perm')
        perm.invoke('all', self._state_update % 'Map reads to %s' % map_to.capitalize())

        # Input files
        hash_v = self._get_index_hash(self._read_length, 'F%d' % self._seed)
        index = File('h%d_%s_F%d_%d.index' % (hash_v, map_to, self._seed, self._read_length))
        reads_txt = File('%s_%s_reads.txt' % (tag, map_to.lower()))

        for i in self._range():
            # Input files
            reads_i = File(reads % i)

            # Output files
            file_type = 'sam' if output_sam else 'mapping'
            path, file_name, ext = GTFAR._get_filename_parts(reads_i.name)
            sam_mapping = '%s_B_%d_%d_%s.%s' % (map_to.upper(), self._seed, self._mismatches, file_name, file_type)
            fastq_out = File('%s_miss_%s%s' % (file_name, map_to, ext))

            # Uses
            perm.uses(reads_i, link=Link.INPUT)
            perm.uses(fastq_out, link=Link.OUTPUT, transfer=False, register=False)
            perm.uses(sam_mapping, link=Link.OUTPUT, transfer=False, register=False)

        # Output files
        log = File('%s_%s.log' % (tag, map_to.upper()))

        # Arguments
        perm.addArguments(index, reads_txt, '--seed F%d' % self._seed, '-v %d' % self._mismatches, '-B', '--printNM')
        perm.addArguments('-u', '-s', '-T %d' % self._read_length)

        if output_sam:
            perm.addArguments('--noSamHeader', '--outputFormat', 'sam')

        perm.setStdout(log)

        # Uses
        perm.uses(index, link=Link.INPUT)
        perm.uses(reads_txt, link=Link.INPUT)
        perm.uses(log, link=Link.OUTPUT, transfer=False, register=False)

        self.adag.addJob(perm)

    def _compute_clip_seed(self, read_length):
        if read_length >= 100:
            anchor = 35
        elif read_length == 75:
            anchor = 25

        return anchor

    def _clip_and_parse_reads(self, reads, tag):
        anchor = self._compute_clip_seed(self._read_length)

        clip_reads = Job(name='clipR')
        clip_reads.invoke('all', self._state_update % 'Generate new splice candidates')

        seed = 'F1'
        mismatches = 1

        # Input files
        hash_v = self._get_index_hash(self._read_length, seed)
        index = File('h%d_%s_F%d_%d.index' % (hash_v, 'GENE', 1, self._read_length))
        reads_txt = File('%s_%s_reads.txt' % (tag, 'gene'))

        for i in self._range():
            # Input files
            reads_i = File(reads % i)

            # Output files
            """
            file_type = 'sam' if output_sam else 'mapping'
            path, file_name, ext = GTFAR._get_filename_parts(reads_i.name)
            sam_mapping = '%s_B_%d_%d_%s.%s' % (map_to.upper(), self._seed, self._mismatches, file_name, file_type)
            fastq_out = File('%s_miss_%s%s' % (file_name, map_to, ext))
            """
            # Uses
            clip_reads.uses(reads_i, link=Link.INPUT)
            #clip_reads.uses(fastq_out, link=Link.OUTPUT, transfer=False, register=False)
            #clip_reads.uses(sam_mapping, link=Link.OUTPUT, transfer=False, register=False)

        # Output files
        log = File('%s_%s.log' % (tag, 'GENE'))

        # Arguments
        clip_reads.addArguments(index, reads_txt, '--seed %s' % seed, '--anchorL %d' % anchor, '-e', '-v %d' % mismatches)
        clip_reads.addArguments('-s', '-u', '--noSamHeader', '--ignoreDummyR %d' % 40, '--ignoreRepeatR %d' % 15)

        clip_reads.setStdout(log)

        # Uses
        clip_reads.uses(index, link=Link.INPUT)
        clip_reads.uses(reads_txt, link=Link.INPUT)
        # TODO: Set Transfer to False for log File.
        clip_reads.uses(log, link=Link.OUTPUT, transfer=True, register=False)

        self.adag.addJob(clip_reads)

    def _parse_alignment(self, input_file, tag):
        parse_alignment = Job(name='parse_alignment')
        parse_alignment.invoke('all', self._state_update % 'Parse alignment')

        # Input files
        input_file = File(input_file)

        # Output files

        vis = File('%s.vis' % input_file.name)

        self._vis_files.append(vis.name)

        # Arguments
        parse_alignment.addArguments(input_file, '--strandRule', self._strand_rule, '--tag', tag)
        parse_alignment.setStdout(vis)

        # Uses
        parse_alignment.uses(input_file, link=Link.INPUT)
        parse_alignment.uses(vis, link=Link.OUTPUT, transfer=False, register=False)

        self.adag.addJob(parse_alignment)


class AnalyzeMixin(object):
    def analyze(self):
        self._analyze()

    def _analyze(self):
        analyze = Job(name='analyze_samfile')
        analyze.invoke('all', self._state_update % 'Analyzing SAM file')

        # Input files
        sam_file = File('%s.sam' % self._prefix)

        # Output files
        genes_counts = File('%s.genes.cnts' % self._prefix)
        features_counts = File('%s.features.cnts' % self._prefix)
        multi_genes_counts = File('%s.multiGenes.cnts' % self._prefix)
        summary_out = File('%s.summary.out' % self._prefix)

        # Arguments
        analyze.addArguments(sam_file, '--prefix', self._prefix)

        # Uses
        analyze.uses(sam_file, link=Link.INPUT)
        analyze.uses(genes_counts, link=Link.OUTPUT, transfer=True, register=False)
        analyze.uses(features_counts, link=Link.OUTPUT, transfer=True, register=False)
        analyze.uses(multi_genes_counts, link=Link.OUTPUT, transfer=True, register=False)
        analyze.uses(summary_out, link=Link.OUTPUT, transfer=True, register=False)

        self.adag.addJob(analyze)


class GTFAR(AnnotateMixin, FilterMixin, IterativeMapMixin, AnalyzeMixin):
    def __init__(self, gtf, genome, prefix, reads, base_dir, bin_dir, read_length=100, mismatches=3,
                 is_trim_unmapped=False, is_map_filtered=False, splice=True, clip_reads=False, strand_rule='Unstranded',
                 dax=None, url=None, email=None, splits=1, adag=None):

        # Reference
        self._gtf = gtf
        self._genome = genome

        # Input
        self._prefix = prefix
        self._reads = reads
        self._read_length = read_length

        if self._read_length > 101:
            self._trims = [50, 75, 100]
        elif self._read_length > 76:
            self._trims = [50, 75]
        elif self._read_length > 51:
            self._trims = [50]

        # Options
        self._mismatches = mismatches
        self._is_trim_unmapped = is_trim_unmapped
        self._is_map_filtered = is_map_filtered
        self._splice = splice
        self._clip_reads = clip_reads
        self._strand_rule = strand_rule

        # Pegasus
        self._base_dir = base_dir
        self._bin_dir = os.path.abspath(bin_dir)
        self._dax = sys.stdout if dax is None else '%s.dax' % dax
        self._url = url
        self._email = email
        self._splits = splits

        def get_range():
            return range(self._splits)

        self._range = get_range

        notifications_conf = os.path.join(self._base_dir, 'config', 'notifications.conf')

        self.adag = AutoADAG('gtfar_%s' % self._prefix) if adag is None else adag

        self._state_update = '%s/gtfar-state-update --id %s --config %s' % (self._bin_dir,
                                                                            self._prefix,
                                                                            notifications_conf)

        self.adag.invoke('all', self._state_update)

        if self._email:
            email_script = '%s/gtfar-email --id %r --from %r --to %r --subject %r --url %r --config %r' % (
                self._bin_dir,
                self._prefix,
                'pegasus-gtfar@localhost.com',
                self._email,
                'GTFAR Workflow finished running',
                self._url,
                notifications_conf)

            self.adag.invoke('at_end', email_script)

        self._state_update = '%s --log-message %%r' % self._state_update

        self._seed = None
        self._vis_files = []

    def validate(self):
        errors = {}

        if self._prefix is None:
            errors['prefix'] = 'Prefix is required'

        self._prefix = str(self._prefix)

        if self._reads is None:
            errors['reads'] = 'Reads needed for options filter'

        elif len(self._reads) == 0:
            errors['reads'] = 'Filter module uses a single read file'

        try:
            self._read_length = int(self._read_length)

            if not 36 <= self._read_length <= 128:
                errors['length'] = 'Read length must be between 50 and 128 (inclusive)'

        except TypeError:
            errors['length'] = 'Read length cannot be None'
        except ValueError:
            errors['length'] = 'Read length must be a valid integer'

        try:
            self._mismatches = int(self._mismatches)

            if not 0 <= self._mismatches <= 8:
                errors['mismatches'] = 'Mismatches must be between 0 and 8 (inclusive)'

        except TypeError:
            errors['mismatches'] = 'Mismatches cannot be None'
        except ValueError:
            errors['mismatches'] = 'Mismatches must be a valid integer'

        if self._clip_reads and self._read_length != 75 and not self._read_length >= 100:
            errors['common'] = 'With clip-reads, read-length must be 75 or >= 100'

        return True if len(errors) == 0 else errors

    def count_and_split_reads(self):
        path, file_name, ext = self.__get_filename_parts(self._reads)
        ext = ext.lower()

        if ext == '.gz':
            pass
        elif ext == '.fastq' or ext == '.fq':
            pass
        else:
            # Invalid extension
            pass

    def _get_index_hash(self, read_length, seed=None, exclude_genome=False):
        hash_k = '%s-%s-%d-%s'
        t = (self._gtf, self._genome, read_length, seed)

        if exclude_genome:
            hash_k = '%s-%d'
            t = (self._gtf, read_length)

        hash_k = hash_k % t
        return mmh3.hash(hash_k, read_length)

    @staticmethod
    def _get_filename_parts(filename):
        path = os.path.dirname(filename)
        file_name, ext = os.path.splitext(filename)
        file_name = os.path.basename(file_name)

        return path, file_name, ext

    def _write_reads_file(self, files_pattern, reads_txt):
        reads_txt = os.path.join(self._base_dir, 'input', reads_txt)
        if not os.path.isdir(os.path.dirname(reads_txt)):
            os.makedirs(os.path.dirname(reads_txt))

        with open(reads_txt, 'w+') as reads:
            for i in self._range():
                reads.write(files_pattern % i + '\n')

    def write_dax(self):
        # Write out reads file
        self._write_reads_file('reads%d_full.fastq', 'full_features_reads.txt')
        self._write_reads_file('reads%d_full_miss_FEATURES.fastq', 'full_genome_reads.txt')
        self._write_reads_file('reads%d_full_miss_FEATURES_miss_GENOME.fastq', 'full_splices_reads.txt')

        if self._clip_reads:
            self._write_reads_file('reads%d_full_miss_FEATURES_miss_GENOME_miss_SPLICES.fastq', 'full_gene_reads.txt')

        if self._is_map_filtered:
            for trim in self._trims:
                self._write_reads_file('reads%%d_%d.fastq' % trim, 'filter%d_features_reads.txt' % trim)
                self._write_reads_file('reads%%d_%d_miss_FEATURES.fastq' % trim, 'filter%d_genome_reads.txt' % trim)
                self._write_reads_file('reads%%d_%d_miss_FEATURES_miss_GENOME.fastq' % trim,
                                       'filter%d_splices_reads.txt' % trim)

        if isinstance(self._dax, str) or isinstance(self._dax, unicode):
            with open(self._dax, 'w') as dax:
                self.adag.writeXML(dax)
        else:
            self.adag.writeXML(sys.stdout)
