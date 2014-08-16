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

from Pegasus.DAX3 import ADAG, Dependency, Job, File, Link, Executable, PFN


class UNIXUtils(object):
    @staticmethod
    def cat(inputs, output):
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
        cat.uses(output, link=Link.OUTPUT, transfer=False, register=False)

        return cat


class AnnotateMixin(object):
    def annotate(self):
        self._annotate(self._read_length)

        if self._is_map_filtered or self._is_trim_unmapped:
            for trim in self._trims:
                self._annotate(trim)

    def _annotate(self, read_length):
        self._annotate_gtf(read_length)

        self._features_index(read_length)
        self._chrs_index(read_length)
        self._splices_index(read_length)
        # self._genes_index(read_length)

    def _annotate_gtf(self, read_length):
        annotate_gtf = Job(name='annotate_gtf')
        annotate_gtf.invoke('all', '%sstate_update.py %r %r %r %r')

        prefix = self._get_index_hash(read_length, exclude_genome=True)

        # Inputs
        gtf = File(self._gtf)
        genome = File(self._genome)

        for i in range(1, 23):
            chr_i = File('chr%d.fa' % i)

            # Uses
            annotate_gtf.uses(chr_i, link=Link.INPUT)

        # Outputs
        features = File('h%s_features.fa' % prefix)
        chrs = File('h%s_chrs.fa' % prefix)
        splices = File('h%s_jxnCands.fa' % prefix)
        genes = File('h%s_geneSeqs.fa' % prefix)

        # Arguments
        annotate_gtf.addArguments(gtf, '-c .', '-p', self._prefix, '-l %d' % read_length)

        # Uses
        annotate_gtf.uses(gtf, link=Link.INPUT)
        annotate_gtf.uses(genome, link=Link.INPUT)
        annotate_gtf.uses(features, link=Link.OUTPUT, transfer=False, register=False)
        annotate_gtf.uses(chrs, link=Link.OUTPUT, transfer=False, register=False)
        annotate_gtf.uses(splices, link=Link.OUTPUT, transfer=False, register=False)
        annotate_gtf.uses(genes, link=Link.OUTPUT, transfer=False, register=False)

        self.adag.addJob(annotate_gtf)

    def _features_index(self, read_length, read_format='fastq', seed='F2'):
        return self._perm_index('features', read_length, read_format=read_format, seed=seed)

    def _chrs_index(self, read_length, read_format='fastq', seed='F2'):
        return self._perm_index('chrs', read_length, read_format=read_format, seed=seed)

    def _splices_index(self, read_length, read_format='fastq', seed='F2'):
        return self._perm_index('jxnCands', read_length, read_format=read_format, seed=seed)

    def _genes_index(self, read_length, read_format='fastq', seed='F1'):
        return self._perm_index('geneSeqs', read_length, read_format=read_format, seed=seed)

    def _perm_index(self, index_type, read_length, read_format='fastq', seed='F2'):
        perm_index = Job(name='perm')
        perm_index.invoke('all', '%sstate_update.py %r %r %r %r')

        prefix = self._get_index_hash(read_length, exclude_genome=True)

        # Input files
        fa_input = File('h%s_%s.fa' % (prefix, index_type))

        # Output files
        hash_v = self._get_index_hash(read_length)
        index = File('h%d_%s_%s_%s.index' % (hash_v, index_type, seed, read_length))

        # Arguments
        perm_index.addArguments(fa_input, '%d' % read_length, '--readFormat', read_format, '--seed', seed)
        perm_index.addArguments('-s', index)

        # Uses
        perm_index.uses(fa_input, link=Link.INPUT)
        # Save this file
        perm_index.uses(index, link=Link.OUTPUT, transfer=False, register=False)

        self.adag.addJob(perm_index)

        return perm_index


class FilterMixin(object):
    def option_filter(self):
        self._option_filter()

    def _option_filter(self):
        option_filter = Job(name='option_filter')
        option_filter.invoke('all', '%sstate_update.py %r %r %r %r')

        # Inputs
        reads = File(self._reads)

        # Outputs
        rejects = File('%s.reject.fastq' % self._prefix)
        adaptor_stats = File('%s.adaptor.stats' % self._prefix)

        # Arguments
        option_filter.addArguments(self._prefix, reads, '%d' % self._read_length)

        # Uses
        option_filter.uses(reads, link=Link.INPUT)

        for i in self._range():
            reads_i = File('reads%d_full.fastq' % i)
            rejects_i = File('reads%d_reject.fastq' % i)
            adaptor_stats_i = File('reads%d.stats' % i)

            for t in self._trims:
                reads_i_t = File('reads%d_%d.fastq' % (i, t))
                option_filter.uses(reads_i_t, link=Link.OUTPUT, transfer=False, register=False)

            option_filter.uses(reads_i, link=Link.OUTPUT, transfer=False, register=False)
            option_filter.uses(rejects_i, link=Link.OUTPUT, transfer=False, register=False)
            option_filter.uses(adaptor_stats_i, link=Link.OUTPUT, transfer=False, register=False)

        option_filter.uses(rejects, link=Link.OUTPUT, transfer=False, register=False)
        option_filter.uses(adaptor_stats, link=Link.OUTPUT, transfer=False, register=False)

        self.adag.addJob(option_filter)


class IterativeMapMixin(object):
    def iterative_map(self):
        self._map_and_parse_reads('reads%d_full.fastq', 'full')
        unmapped = 'reads%s_full_unmapped.fastq'

        for trim in self._trims:
            self._read_length = trim

            if (self._is_trim_unmapped):
                self._map_and_parse_reads(unmapped, 'trim%d' % trim)
                unmapped = 'reads*unmapped.fastq'

            if (self._is_map_filtered):
                self._map_and_parse_reads('reads%%d_%d.fastq' % trim, 'filter%d' % trim)

        cat = UNIXUtils.cat(self._vis_files, '%s.sam' % self._prefix)
        cat.invoke('all', '%sstate_update.py %r %r %r %r')
        self.adag.addJob(cat)

        if self._is_trim_unmapped:
            cat = UNIXUtils.cat([2], '%s.unmapped.fastq')
        else:
            cat = UNIXUtils.cat(['reads%d_full_unmapped.fastq' % i for i in self._range()], '%s.unmapped.fastq')

        cat.invoke('all', '%sstate_update.py %r %r %r %r')
        self.adag.addJob(cat)

    def _map_and_parse_reads(self, reads, tag):
        self._setup_perm_seeds()

        self._map_and_parse_reads_to_features(reads, tag)

        path, file_name, ext = GTFAR._get_filename_parts(reads)
        reads = '%s_miss_FEATURES%s' % (file_name, ext)
        self._map_and_parse_reads_to_genome(reads, tag)

        if self._splice:
            reads = '%s_miss_GENOME%s' % os.path.splitext(reads)
            self._map_and_parse_reads_to_splices(reads, tag)

    def _setup_perm_seeds(self):
        self._seed = self._mismatches

        if self._read_length > 64:
            self._seed = (self._mismatches + 1) / 2

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
        perm.invoke('all', '%sstate_update.py %r %r %r %r')

        # Input files
        hash_v = self._get_index_hash(self._read_length)
        index = File('h%d_%s_F%d_%d.index' % (hash_v, index_type, self._seed, self._read_length))
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
        perm.addArguments('-u', '-s', '-T %d' % self._read_length, '--log', log)

        if output_sam:
            perm.addArguments('--noSamHeader', '--outputFormat', 'sam')

        # Uses
        perm.uses(index, link=Link.INPUT)
        perm.uses(reads_txt, link=Link.INPUT)
        perm.uses(log, link=Link.OUTPUT, transfer=False, register=False)

        self.adag.addJob(perm)

    def _parse_alignment(self, input_file, tag):
        parse_alignment = Job(name='parse_alignment')
        parse_alignment.invoke('all', '%sstate_update.py %r %r %r %r')

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


class GTFAR(AnnotateMixin, FilterMixin, IterativeMapMixin):
    def __init__(self, gtf, genome, prefix, reads, base_dir, read_length=100, mismatches=3, is_trim_unmapped=False,
                 is_map_filtered=False, splice=False, strand_rule='Unstranded', dax=None, url=None, email=None,
                 adag=None):

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
        self._splice = True
        self._strand_rule = strand_rule

        # Pegasus
        self._base_dir = base_dir
        self._dax = sys.stdout if dax is None else '%s.dax' % (dax)
        self._url = url
        self._email = email

        def get_range():
            return range(0, 1)

        self._range = get_range

        self.adag = ADAG('gtfar_%s' % self._prefix) if adag is None else adag

        self._state_update = '%sstate_update.py %%r %%r %%r %%r' % '/bin/dir'
        self._email_script = '%sgtfar-email %%r %%r %%r %%r' % '/bin/dir'

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

    def _get_index_hash(self, read_length, exclude_genome=False):
        hash_k = '%s-%s-%d'
        t = (self._gtf, self._genome, read_length)

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
        self._write_reads_file('reads%d_full.fastq', 'full_feature_reads.txt')
        self._write_reads_file('reads%d_full_miss_FEATURES.fastq', 'full_genome_reads.txt')
        self._write_reads_file('reads%d_full_miss_FEATURES_miss_GENOME.fastq', 'full_splices_reads.txt')

        if self._is_map_filtered:
            for trim in self._trims:
                self._write_reads_file('reads%%d_%d.fastq' % trim, 'filter%d_feature_reads.txt' % trim)
                self._write_reads_file('reads%%d_%d_miss_FEATURES.fastq' % trim, 'filter%d_genome_reads.txt' % trim)
                self._write_reads_file('reads%%d_%d_miss_FEATURES_miss_GENOME.fastq' % trim,
                                       'filter%d_splices_reads.txt' % trim)

        if isinstance(self._dax, str):
            with open(self._dax, 'w') as dax:
                self.adag.writeXML(dax)
        else:
            self.adag.writeXML(sys.stdout)
