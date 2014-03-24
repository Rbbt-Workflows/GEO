$LOAD_PATH.unshift(File.join(File.dirname(__FILE__), '../..', 'lib'))
$LOAD_PATH.unshift(File.dirname(__FILE__))
require 'test/unit'
require 'expression/matrix'

class TestClass < Test::Unit::TestCase
  def test_true
    assert true
  end
end

