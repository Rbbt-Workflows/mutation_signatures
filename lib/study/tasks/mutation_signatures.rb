require 'sample/tasks/mutation_signatures' 

if Module === Study and Workflow === Study
  module Study

    dep Sample, :mutation_signature, :compute => :bootstrap do |jobname,options|
      study = Study.setup(jobname.dup)
      study.genotyped_samples.collect do |sample|
        sample.mutation_signature(:job)
      end
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

    tasks[:nmf_features] = MutationSignatures.tasks[:nmf_features]
    tasks[:nmf_profile_plots] = MutationSignatures.tasks[:nmf_profile_plots]
    task_dependencies[:nmf_features] = MutationSignatures.task_dependencies[:nmf_features]
    task_dependencies[:nmf_profile_plots] = MutationSignatures.task_dependencies[:nmf_profile_plots]


    tasks[:em_features] = MutationSignatures.tasks[:em_features]
    tasks[:em_profile_plots] = MutationSignatures.tasks[:em_profile_plots]
    task_dependencies[:em_features] = MutationSignatures.task_dependencies[:em_features]
    task_dependencies[:em_profile_plots] = MutationSignatures.task_dependencies[:em_profile_plots]
  end
end
