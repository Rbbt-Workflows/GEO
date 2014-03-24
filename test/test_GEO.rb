$LOAD_PATH.unshift(File.join(File.dirname(__FILE__), '..', 'lib'))
$LOAD_PATH.unshift(File.dirname(__FILE__))
require 'test/unit'
require 'GEO'

class TestGEO < Test::Unit::TestCase
  def test_query
    assert_equal ["GDS1761"], GEO.query("NCI60[Title]")
  end
end

