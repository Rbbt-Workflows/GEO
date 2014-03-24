class Matrix
  def differential(main, contrast, path)
    if Array === main and Array === contrast
      main_samples, contrast_samples = main, contrast
    else
      main_samples, contrast_samples = comparison labels, main, contrast
    end

    log2 = value_type == "count"
    GE.analyze(data_file, main_samples, contrast_samples, log2, path, format)
  end
end
