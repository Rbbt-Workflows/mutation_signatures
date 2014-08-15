require 'sample/tasks/mutation_signatures' 

module Study

  dep Sample, :mutation_signature do |jobname,options|
    Study.samples(jobname).collect do |sample|
      next unless sample.has_genotype?
      sample = [jobname, sample] * ":"
      Sample.job(:mutation_signature, sample)
    end.compact
  end
  task :cohort_signatures => :tsv do

    samples = dependencies.collect{|dep| dep.clean_name }
    fields = samples
    signatures = TSV.setup({}, :key_field => "Change", :fields => fields, :type => :list, :cast => :to_i)
    Step.wait_for_jobs dependencies
    dependencies.each do |dep|
      sample = dep.clean_name
      signature = dep.load
      sample_pos = samples.index sample
      signature.each do |change,count|
        signatures[change] ||= [0] * samples.length
        signatures[change][sample_pos] = count
      end
    end

    signatures
  end
end
