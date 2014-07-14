require 'rbbt'
require 'rbbt/workflow'
require 'rbbt/util/misc'

Workflow.require_workflow "Genomics"

require 'rbbt/sources/organism'
require 'rbbt/entity/genomic_mutation'
require 'set'
require 'bio'
require File.join(File.dirname(__FILE__), 'exome_positions')

module MutationSignatures
  input :mutations, :array, "Mutations"
  input :organism, :string, "Organism code"
  def self.context(mutations, organism)
    chr_mutations = {}
    mutations.each do |mutation|
      chr, *rest = mutation.split(":")
      chr_mutations[chr] ||= []
      chr_mutations[chr] << mutation
    end

    result = TSV.setup({}, :key_field => "Genomic Mutation", :fields => ["Context change"], :type => :single)
    chr_mutations.each do |chr, list|
      begin
        chr = chr.sub('chr','')
        Organism[organism]["chromosome_" << chr].open do |chr_file|
          positions = list.collect{|mutation| chr, pos, *rest = mutation.split(":"); pos.to_i}
          position_context = {}
          positions.sort.each do |pos|
            chr_file.seek(pos-2)
            context = chr_file.read(3)
            position_context[pos] = context
          end
          list.each do |mutation|
            chr, pos, *rest =mutation.split(":")
            pos = pos.to_i
            context = position_context[pos]
            result[mutation] = context
          end
        end
      rescue Exception
        Log.debug("Skipping #{ chr }: #{$!.message}")
        next
      end
    end
    result
  end
  task :context => :tsv
  export_asynchronous :context

  Rbbt.claim Rbbt.software.opt.EMu, :install, Rbbt.share.install.software.EMu.find

  input :exome_only, :boolean, "Limit to exome bases", false
  task :mutation_oportunities => :tsv do |exome_only|
    organism = name.split("_").first.sub('-', '/')

    contexts = Set.new

    chr_contexts = {}
    Organism.chromosomes(organism).tsv.keys.each do |chromosome|
      begin
        chr_file = Organism[organism]['chromosome_' << chromosome].open
      rescue
        Log.debug("Error processing mutation oportunities: " << $!.message)
        next
      end

      if exome_only
        log(:exome_only, "Collapsing exon ranges: #{ chromosome }")
        exon_ranges = MutationSignatures.exon_ranges_for_chr(organism)
        exon_ranges.filter
        exon_ranges.add_filter("field:Chromosome Name", chromosome)
        collapsed_ranges = Misc.collapse_ranges(exon_ranges.sort_by("Exon Chr Start"){|v,k,all| v.to_i}.collect{|e,range| (range[0].to_i..range[1].to_i)})
        exon_ranges.pop_filter
        exon_ranges.reset_filters


        str = ""
        log(:exome_only, "Subsetting chromosome sequence: #{ chromosome }")
        collapsed_ranges.each do |range|
          chr_file.seek(range.begin)
          str << chr_file.read(range.end - range.begin + 1)
        end
        chr_file = StringIO.new str
      end

      chr_context = {}
      last = [nil,chr_file.getc,chr_file.getc]
      log(:context, "Counting context events: #{ chromosome } (size #{chr_file.size})")
      while c = chr_file.getc
        last.shift
        last.push c
        case last[1]
        when "T", "C"
          context = last * ""
        when "A", "G"
          context = Bio::Sequence::NA.new(last * "").complement.upcase.reverse
        else
          next
        end
        contexts << context
        chr_context[context] ||= 0
        chr_context[context] += 1
      end
      chr_contexts[chromosome] = chr_context
    end

    contexts = contexts.sort
    tsv = TSV.setup({}, :key_field => "Chromosome Name", :fields => contexts, :type => :list, :cast => :to_i)

    chr_contexts.each do |chromosome, values|
      tsv[chromosome] = values.values_at *contexts
    end

    tsv
  end

  input :mutations, :array, "Mutations"
  input :organism, :string, "Organism code"
  input :watson, :boolean, "Mutations reported in Watson strand", true
  task :mutation_context => :tsv do |mutations, organism, watson|

    if not watson
      log :watson, "Turning mutations to watson"
      mutations = Sequence.job(:to_watson, name, :mutations => mutations, :organism => organism).run
    end

    log :mutation_context, "Getting surrounding bases for each mutation"
    mutation_context = MutationSignatures.context(mutations, organism)

    log :changes, "Turning changes into context changes"
    changes = mutations.collect do |m| 
      base = m.split(":")[2]
      next unless %w(A C T G).include? base

      begin
        context = mutation_context[m]
      rescue
        Log.debug("Exception extracting context: " << $!.message)
        next
      end
      next if context.nil? or base == context[1]

      if context[1] == "G" or context[1] == "A"
        context = Bio::Sequence::NA.new(context).complement.upcase.reverse
        base = Bio::Sequence::NA.new(base).complement.upcase
      end

      [context, base] * ">"
    end

    log :counts, "Counting changes"
    counts = Misc.counts(changes.compact).dup

    counts = TSV.setup(counts, :key_field => "Change", :fields => ["Count"], :type => :single, :cast => :to_i)

    counts
  end
  export_asynchronous :mutation_context


  #input :mutations, :array, "Mutations"
  #input :organism, :string, "Organism code"
  #input :watson, :boolean, "Mutations reported in Watson strand", true
  #task :mutation_context => :tsv do |mutations, organism, watson|
  #  GenomicMutation.setup(mutations, name, organism, watson)

  #  changes = mutations.collect do |m| 
  #    base = m.base
  #    next unless %w(A C T G).include? base

  #    begin
  #      context = m.bases_in_range(-1..1)
  #    rescue
  #      Log.debug("Exception extracting context: " << $!.message)
  #      next
  #    end
  #    next if context.nil? or base == context[1]


  #    if context[1] == "G" or context[1] == "A"
  #      context = Bio::Sequence::NA.new(context).complement.upcase.reverse
  #      base = Bio::Sequence::NA.new(base).complement.upcase
  #    end

  #    [context, base] * ">"
  #  end
  #  counts = Misc.counts(changes.compact).dup

  #  counts = TSV.setup(counts, :key_field => "Change", :fields => ["Count"], :type => :single, :cast => :to_i)

  #  counts
  #end
  #export_asynchronous :mutation_context

  input :cohort, :tsv, "Genomic Sample and Genomic Mutations. Example row: 'Sample01{TAB}10:12345678:A{TAB}X:5454322:G'"
  input :organism, :string, "Organism code"
  input :watson, :boolean, "Mutations reported in Watson strand", true
  task :cohort_signatures => :tsv do |cohort, organism, watson|
    samples = cohort.keys.sort

    signatures = TSV.setup({}, :key_field => "Change", :fields => samples, :type => :list, :cast => :to_i)
    cohort.each do |sample, mutations|
      signature = MutationSignatures.job(:mutation_context, nil, :mutations => mutations.flatten, :organism => organism, :watson => watson).exec

      sample_pos = samples.index sample
      signature.each do |change,count|
        signatures[change] ||= [0] * samples.length
        signatures[change][sample_pos] = count
      end
    end

    signatures
  end
  export_asynchronous :cohort_signatures

  dep :cohort_signatures
  dep do |jobname, inputs| MutationSignatures.job(:mutation_oportunities, inputs[:organism].sub('/', '-'), :exome_only => inputs[:exome_only]) end
  input :scale, :select, "Scale", 'none', :select_options => [:none, :sample, :channel, :both]
  input :exome_only, :boolean, "Limit to exome bases", false
  task :signature_plots => :binary do |scale|
    signatures = step(:cohort_signatures).load
    opportunities = step(:mutation_oportunities).load

    height = 5 * signatures.size
    height = 50 if height > 50
    TmpFile.with_file(opportunities.to_s, false) do |opp|
      FileUtils.mkdir files_dir unless File.exists? files_dir
      script =<<-EOF
source('#{Rbbt.share.R["plots.R"].find(:lib)}')
opportunities <- rbbt.tsv('#{opp}');
channel_counts = apply(opportunities, 2, sum)
p <- signature_heatmap(data, channel_counts, #{R.ruby2R(scale)})
ggsave('#{file('heatmap.png')}', p, height=#{height}, limitsize=FALSE)
data  = NULL
      EOF

      signatures.R script
    end

    Open.read(file('heatmap.png'), :mode => "rb")
  end
  export_asynchronous :signature_plots

  #{{{ NMF

  dep :cohort_signatures
  input :k, :integer, "Number of features"
  task :nmf_features => :tsv do |k|
    signatures = step(:cohort_signatures).load

    FileUtils.mkdir files_dir unless File.exists? files_dir
    script =<<-EOF

library(NMF)

nmf = nmf(data, #{ k }, method='lee');
str(data)
names = sapply((1:#{ k }), function(n){paste("Factor",n,sep=" ")})

w = nmf@fit@W
colnames(w) = names
h = nmf@fit@H
rownames(h) = names

rbbt.tsv.write('#{file('factor_composition.tsv')}', w)
rbbt.tsv.write('#{file('sample_factors.tsv')}', h)

data = w
    EOF

    signatures.R script
  end
  export_asynchronous :nmf_features


  dep :nmf_features
  task :nmf_profile_plots => :binary do
    nmf_factors = step(:nmf_features).file("factor_composition.tsv").tsv
    nmf_samples = step(:nmf_features).file("sample_factors.tsv").tsv

    FileUtils.mkdir files_dir unless File.exists? files_dir

    height = nmf_factors.size * 0.1
    height = 50 if height > 50
    width = 7
    script =<<-EOF
source('#{Rbbt.share.R["plots.R"].find(:lib)}')
p <- factor_profile_plot(data)
ggsave('#{file('factor_profile.png')}', p, height=#{height}, width=#{width}, limitsize=FALSE)
data  = NULL
    EOF
    #nmf_factors.R script

    height = nmf_factors.size * 0.1
    height = 50 if height > 50
    width = nmf_samples.fields.length 
    width = 50 if width > 50
    script =<<-EOF
source('#{Rbbt.share.R["plots.R"].find(:lib)}')
p <- sample_profile_plot(data)
ggsave('#{file('sample_profile.png')}', p, height=#{height}, width=#{width},  limitsize=FALSE)
data  = NULL
    EOF
    ovary_samples = nmf_samples.fields.select{|f| f =~ /squamous/i}
    nmf_samples.slice(ovary_samples).R script


    Open.read(file('factor_profile.png'), :mode => "rb")
  end
  export_asynchronous :nmf_profile_plots

  #{{{ EM

  dep :cohort_signatures
  dep do |jobname, inputs| MutationSignatures.job(:mutation_oportunities, inputs[:organism].sub('/', '-'), :exome_only => inputs[:exome_only]) end
  input :k, :integer, "Number of spectra (0 = use BIC)", 0
  input :exome_only, :boolean, "Limit to exome bases", false
  task :em_features => :integer do |k,eo|
    signatures = step(:cohort_signatures).load
    oportunities = step(:mutation_oportunities).load
    
    channels = signatures.keys.select{|c| c =~ /^[ACTG>]*$/}.sort_by do |c|
      context, base = c.split ">"
      [context[1], base, context] * ":"
    end

    channel_contexts  = Misc.process_to_hash(channels){|channels| channels.collect{|c| c.split(">").first}}

    global_oportunities = {}
    channels.collect do |channel|
      global_oportunities[channel] = 0
      oportunities.each do |chr, values|
        context = channel_contexts[channel]
        global_oportunities[channel] += values[context]
      end
    end

    samples = signatures.fields

    FileUtils.mkdir files_dir unless File.exists? files_dir
    TmpFile.with_file(global_oportunities.values_at(*channels) * "\t", false) do |opp|
      TmpFile.with_file(signatures.slice(samples).values_at(*channels).transpose.collect{|l| l * "\t"} * "\n", false) do |mut|
        CMD.cmd("#{Rbbt.software.opt.EMu.EMu.find} --mut #{ mut } --opp #{opp} --pre #{File.join(files_dir, 'result')} #{k.nil? or k.to_i == 0 ? "" : "--force #{ k }"}")
      end
    end

    set_info :channels, channels
    set_info :samples, samples

    Dir.glob(File.join(files_dir, 'result_*_ml_spectra.txt')).first.match(/(\d+)_ml_spectra/)[1].to_i
  end
  export_asynchronous :em_features

  dep :em_features
  task :em_profile_plots => :binary do
    em_features = step(:em_features)
    channels = em_features.info[:channels]
    samples = em_features.info[:samples]

    num_spectra = em_features.load

    assignments = Dir.glob(File.join(em_features.files_dir,"result_*_assigned.txt")).first
    activities  = Dir.glob(File.join(em_features.files_dir,"result_*_map_activities.txt")).first
    spectra     = Dir.glob(File.join(em_features.files_dir,"result_*_ml_spectra.txt")).first


    FileUtils.mkdir files_dir unless File.exists? files_dir
    script =<<-EOF
source('#{Rbbt.share.R["plots.R"].find(:lib)}')

assignments = as.data.frame(t(read.table('#{assignments}')))
activities  = as.data.frame(t(read.table('#{activities}')))
spectra     = as.data.frame(t(read.table('#{spectra}')))

channels = #{R.ruby2R(channels)}
samples = #{R.ruby2R(samples)}

spectra.names = sapply((1:#{num_spectra}), function(n){ paste("Spectra", n, sep=" ")})

rownames(assignments) <- spectra.names
colnames(assignments) <- samples

rownames(activities) <- spectra.names
colnames(activities) <- samples

rownames(spectra) <- channels
colnames(spectra) <- spectra.names

p <- factor_profile_plot(spectra)
ggsave('#{file('factor_profile.png')}', p)

p <- sample_profile_plot(activities)
ggsave('#{file('sample_profile_activities.png')}', p)

p <- sample_profile_plot(assignments)
ggsave('#{file('sample_profile_assignments.png')}', p)
    EOF

    R.run script
 
    Open.read(file('factor_profile.png'), :mode => 'rb')
  end
  export_asynchronous :em_profile_plots

end
