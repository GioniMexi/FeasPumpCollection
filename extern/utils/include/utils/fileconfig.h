/**
 * @file fileconfig.h
 * @brief Config files reader/writer
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 */

#ifndef FILECONFIG_H
#define FILECONFIG_H

#include <map>
#include <string>
#include <memory>
#include "str_utils.h"

namespace dominiqs {


/**
 * Config File Reader/Writer
 */

class FileConfig
{
public:
	/** File I/O */
	//@{
	/**
	 * Load config file
	 * @param fileName file to load
	 * @param merge merges with current contents if set to true, replaces otherwise
	 */
	bool load(const std::string& fileName, bool merge = true);
	/**
	 * Save configuration to file
	 * @param fileName file to save
	 */
	void save(const std::string& fileName) const;
	//@}

	/** typed wrapper around get/set */
	//@{
	template<typename T> T get(const std::string& entry, const T& defValue) const
	{
		return from_string<T>(getAsString(entry, std::to_string(defValue)));
	}
	template<typename T> void set(const std::string& entry, const T& value)
	{
		setAsString(entry, std::to_string(value));
	}
	//@}
private:
	using EntryMap = std::map<std::string, std::string>;
	EntryMap entries;
	/** std::string get/set */
	//@{
	/**
	 * Get the value of a parameter
	 * @param @entry parameter name
	 * @param @defValue default value if parameter not found
	 * @return the parameter value (or @param defValue if not found)
	 */
	std::string getAsString(const std::string& entry, const std::string& defValue) const;
	/**
	 * Set the value of a parameter. Will create one if it does not exists yet
	 * @param @entry parameter name
	 * @param @value parameter value
	 */
	void setAsString(const std::string& entry, const std::string& value);
	//@}
};

/**
 * template specialization of get/set for std::string
 */

template<> std::string FileConfig::get(const std::string& entry, const std::string& defValue) const;

template<> void FileConfig::set(const std::string& entry, const std::string& value);

/**
 * Global Config File
 */

FileConfig& gConfig();

} // namespace dominiqs

#endif /* FILECONFIG_H */
